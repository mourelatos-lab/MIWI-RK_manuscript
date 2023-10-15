#!/usr/bin/env Rscript
#
# Get differential expression of piRNA sequences "annotated" by their seed (2-7 nt)
# Plus, take a list of "special" piRNA-gene (for example those annotated to target genes in GTBuster = piRNA-gene pairs) and
#   add this as an additional requirement
# The piRNAs we have at the end as an expression "unit" are piRNAs with the same seed and targeting the same gene
# We can have piRNAs with the same seed but targeting a different gene or no gene at all. These piRNAs are put in
#   a separate category
#
################################################################################

suppressPackageStartupMessages(library("optparse"))
library("stringr")
library("reshape2") # Used only a couple of times
suppressPackageStartupMessages(library("GenomicRanges"))
library("RColorBrewer")
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("edgeR"))

################################################################################

read_bed <- function(bed) {
  tab <- read.table(bed, header = F, sep = "\t", stringsAsFactors = F)
  colnames(tab)[1:6] <- c("chr", "start", "end", "name", "score", "strand")
  return(tab)
}

coords_to_name <- function(bed) {
  bed$name <- paste(bed$name, bed$chr, bed$start, bed$end, bed$strand, sep = ",")
  return(bed)
}

bed_to_granges <- function(bed) {
  annotGr <- makeGRangesFromDataFrame(bed,
    seqnames.field = "chr", start.field = "start", end.field = "end", strand.field = "strand",
    keep.extra.columns = T,
    starts.in.df.are.0based = T, ignore.strand = F,
    seqinfo = Seqinfo(seqnames = as.character(unique(bed$chr)))
  )
  return(annotGr)
}

# Merge duplicated vector to unique, "|" separated vector
merge_vector <- function(x) {
  unique(x) %>% paste(collapse = "|")
}

################################################################################

option_list <- list(
  make_option(
    c("-i", "--ifile"),
    type = "character",
    help = "Input bed/sequence tsv for all samples.", metavar = "File"
  ),
  make_option(
    c("-d", "--design"),
    type = "character", default = NULL,
    help = "Description file for differential expression.", metavar = "File"
  ),
  make_option(
    c("-o", "--odir"),
    type = "character",
    help = "Output directory where all tables and plots will be saved.", metavar = "Directory"
  ),
  make_option(
    c("-c", "--comparcond1"),
    type = "character",
    help = "Condition 1 to compare."
  ),
  make_option(
    c("-m", "--comparcond2"),
    type = "character",
    help = "Condition 2 to compare."
  ),
  make_option(
    c("-r", "--refcond"),
    type = "character", default = NULL,
    help = "Reference condition. Default: use --comparcond1."
  ),
  make_option(
    c("-s", "--samples"),
    type = "character",
    help = "List of samples to compare separated by a single comma."
  ),
  make_option(
    c("-p", "--pval"),
    type = "double", default = 0.001,
    help = "Adjusted p-value to use as cutoff. Default: 0.001."
  ),
  make_option(
    c("-f", "--fc"),
    type = "double", default = 2.0,
    help = "Fold-change to use as cutoff. Default: 2."
  ),
  make_option(
    c("-l", "--pirna_list"),
    type = "character", default = NULL,
    help = "List of special piRNA-gene pairs for making sub-seed families (like those from GTBuster). Default: NULL."
  ),
  make_option(
    c("-t", "--type"),
    type = "character", default = "seq",
    help = "Summarized counts by positions or by sequence [pos|seq]. Default: seq"
  ),
  make_option(
    c("-b", "--byfeature"),
    type = "logical", default = "FALSE", action = "store_true",
    help = "Calculate differential expression by mapped feature instead of position/sequence. Default: FALSE", metavar = "Boolean"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))
print(opt)

intab <- opt$ifile
design <- opt$design
odir <- opt$odir
condsToCompare <- c(opt$comparcond1, opt$comparcond2)
if (is.null(opt$refcond)) {
  refcond <- opt$comparcond1
} else {
  refcond <- opt$refcond
}
samples <- strsplit(opt$samples, split = ",")[[1]]
pval <- opt$pval
fc <- log2(opt$fc)
top <- opt$top
type <- opt$type
byfeature <- opt$byfeature
pirna_list <- opt$pirna_list
seed_reg <- c(2, 7) # start and end of the seed region; hardcoded for now

################################################################################

dir.create(odir, recursive = T)

coldata <- read.table(file = design, header = T, sep = "\t", stringsAsFactors = F)

coldata <- coldata %>%
  filter(sample %in% samples)
coldata$condition <- factor(coldata$condition, levels = unique(coldata$condition))
coldata$batch <- factor(coldata$batch, levels = unique(coldata$batch))
coldata$condition <- relevel(coldata$condition, ref = refcond)
coldata$batch <- relevel(coldata$batch, ref = levels(coldata$batch)[1])

tab_main <- rio::import(file = intab, format = "tsv")

# Get only selected samples and their coords
tab_main <- tab_main[
  ,
  sort(unlist(lapply(c("name", samples), function(x) grep(x, colnames(tab_main)))))
]
tab_main[, samples][is.na(tab_main[, samples])] <- 0 # Replace NA in expression with 0
tab_main <- tab_main[rowSums(tab_main[, samples] == 0) < length(samples), ] # Get only rows with expression in at least one sample
tab_main[, paste0(samples, "_coords")][is.na(tab_main[, paste0(samples, "_coords")])] <- "." # Replace NA in coords with .

# Load selected piRNAs
if (!is.null(pirna_list)) {
  pirnas <- read.table(pirna_list, header = F, stringsAsFactors = F)
  colnames(pirnas) <- c("name", "gene_id")

  tab_main <- tab_main %>%
    left_join(pirnas)

  # Add whether these piRNA are selected/not selected by GTBusters to be paired with gene
  #   If we do simple seed-counting we might be counting piRNAs which are not predicted to
  #   be interacting with a gene but still have the same seed
  # This will help us to keep piRNA counts even if they cannot be associated with any gene
  tab_main$gene_id[is.na(tab_main$gene_id)] <- "notselect"

  control_cnt <- colSums(tab_main[, 2:5]) # Temp. save control count to check if we didn't lose anything
}

# Extract seed
tab_main$pirna_seed <- substr(tab_main$name, min(seed_reg), max(seed_reg))

# Summarize by seed instead of sequence - separate "selected" piRNAs from non-selected

tab_main <- tab_main %>%
  filter(!is.na(gene_id)) %>%
  group_by(gene_id, pirna_seed) %>%
  mutate(
    name_merged = merge_vector(name),
    across(ends_with("_coords"),
      .fns = list("merged" = ~ merge_vector(.)),
      .names = "{col}_{fn}"
    )
  ) %>%
  ungroup()

# Save coords and piRNA-piRNA seed links
coords <- tab_main %>%
  select(gene_id, name, name_merged, pirna_seed, contains("_coords")) %>%
  distinct()

tab_main <- tab_main %>%
  select(name_merged, !contains("_coords")) %>%
  select(-name, -gene_id, -pirna_seed) %>%
  group_by(name_merged) %>%
  summarise(across(all_of(samples), ~ sum(.x, na.rm = TRUE)))

### Some control check counts
# How many merged "piRNA names" do we have? Should be 0
dupl_merged_pirnas <- table(tab_main$name_merged)
if (sum(dupl_merged_pirnas > 1) != 0) {
  stop("We have problem with duplicated merged piRNA names, this should not happen! Please check.")
}

# Check if all the merged names are in coordinates file
if (sum(!(coords$name_merged %in% tab_main$name_merged)) != 0) {
  stop("We have problem with some merged piRNA names not being present in \"final\" de count input file. This should not happen! Please check.")
}

colSums(tab_main[, 2:5]) == control_cnt # Compare control counts so we see we didn't lose anything. Should be all TRUE

# Get suffix for output files
suffix <- "-bySeq"

# Reorder count file
tab_main <- tab_main %>%
  rename(name = name_merged) %>%
  select(name, !!samples) %>%
  as.data.frame()

rownames(tab_main) <- tab_main$name
tab_main$name <- NULL

tab_main <- tab_main %>%
  select(!!coldata$sample)

# Make sure coldata and counts table have the same order
if (sum((colnames(tab_main) == coldata$sample) == F) != 0) {
  stop("Count table and coldata table don't have the same order, please check and rerun")
} else { # rename to a nice name
  colnames(tab_main) <- coldata$name
}

# Run DE for each group separately
annot_group <- "all"
print(annot_group)

tab <- tab_main

### "Differential expression"
d <- DGEList(counts = tab, group = coldata$condition) # edgeR DGE object

design <- model.matrix(~ 0 + batch + condition, data = coldata)

d <- estimateDisp(d, design)

fit <- glmQLFit(d, design, dispersion = d$tagwise.dispersion)

if (!length(grep(paste0("condition", condsToCompare[1], "\\b"), colnames(fit$design)))) { # If I cannot find condsToCompare[1] (usually intercept) set contrast as 1 for condsToCompare2
  my.contrasts <- makeContrasts( # Should create contrasts; if we have intercept and we want to compare coef=2 it should be the same as contrast=c(0, -1, 0) but it might be misunderstood because there is actually no contrast https://www.biostars.org/p/102036/
    postvspre = paste0(colnames(fit$design)[grep(paste0("condition", condsToCompare[2], "\\b"), colnames(fit$design))]), levels = fit$design
  ) # https://stackoverflow.com/questions/26813667/how-to-use-grep-to-find-exact-match
} else { # Else make proper contrast
  my.contrasts <- makeContrasts( # Should create contrasts; if we have intercept and we want to compare coef=2 it should be the same as contrast=c(0, -1, 0) but it might be misunderstood because there is actually no contrast https://www.biostars.org/p/102036/
    postvspre = paste0(colnames(fit$design)[grep(paste0("condition", condsToCompare[2], "\\b"), colnames(fit$design))], "-", colnames(fit$design)[grep(paste0("condition", condsToCompare[1], "\\b"), colnames(fit$design))]), levels = fit$design
  ) # https://stackoverflow.com/questions/26813667/how-to-use-grep-to-find-exact-match
}
colnames(my.contrasts) <- "contrast"
lrt <- glmQLFTest(fit, contrast = my.contrasts[, "contrast"]) # If we have 3 conditions and want to compare 3 vs 1 we set contrast=c(0, -1, 1), if we want to compare 3 vs 1 or 2 vs 1 we just set coef=3 or coef=2, respectively; some more examples of contrast https://www.biostars.org/p/110861/
lrt_tab <- topTags(lrt, n = nrow(lrt$table), adjust.method = "BH", sort.by = "p.value")$table # Extract all genes

cond_colours <- brewer.pal(length(unique(d$samples$group)), "Paired")[d$samples$group]
names(cond_colours) <- d$samples$group

# With batch effect removal
A <- aveLogCPM(d)
d2 <- d[A > 1, ]
d2 <- calcNormFactors(d2)
logCPM <- cpm(d2, log = TRUE, prior.count = 0.1)
logCPMc <- limma::removeBatchEffect(logCPM, coldata$batch)

pdf(file = paste0(odir, "/MDS_plot_batchEffect", ".", annot_group, suffix, ".pdf"), width = 10)
  par(mfrow = c(1, 2), oma = c(1, 0, 0, 0), xpd = NA)
  plotMDS(logCPM, col = cond_colours, main = "MDS without sample pairing (logCPM)", font = 2)
  plotMDS(logCPMc, col = cond_colours, main = "MDS with sample pairing (logCPM)", font = 2)
  legend(-1, -0.45, levels(d$samples$group), fill = cond_colours[levels(d$samples$group)], cex = 0.6)
dev.off()

de.genes <- rownames(lrt_tab)[(lrt_tab$FDR < pval) & (abs(lrt_tab$logFC) >= fc)]

png(file = paste0(odir, "/edgeR_MAplot_", condsToCompare[1], "_vs_", condsToCompare[2], ".", annot_group, suffix, ".png"), height = 960, width = 960)
  plotSmear(lrt, de.tags = de.genes, main = "FC Plot With Tagwise Dispersion")
  abline(h = c(-fc, fc), col = "dodgerblue")
dev.off()

lrt_tab$name <- rownames(lrt_tab)
colnames(tab) <- paste0(colnames(tab), "_rawCounts")
tab$name <- rownames(tab)

normcounts <- as.data.frame(cpm(d))
colnames(normcounts) <- paste0(colnames(normcounts), "_normCounts")
normcounts$name <- rownames(normcounts)

lrt_tab <- lrt_tab %>%
  left_join(normcounts) %>%
  left_join(tab)

# Add expression bins
# Make 20 even groups by expression
quants <- quantile(lrt_tab$logCPM, probs = seq(0, 1, 0.05)) # just print the quartiles
quants
if (any(duplicated(quants))) { # If we have any duplicated quants just add a tiny number
  print("There are duplicated quantile limits, adding a tiny number.")
  quants[duplicated(quants)] <- quants[duplicated(quants)] + runif(sum(duplicated(quants)), 1e-9, 1e-8)
}
lrt_tab <- within(lrt_tab, AbundanceBin <- as.integer(cut(lrt_tab$logCPM, quants, include.lowest = TRUE), include.lowest = TRUE)) # https://stackoverflow.com/questions/7508229/how-to-create-a-column-with-a-quartile-rank

lrt_tab <- lrt_tab %>%
  mutate(change = case_when(
    (FDR < !!pval) & (logFC >= !!fc) ~ "Up",
    (FDR < !!pval) & (logFC <= (-!!fc)) ~ "Down",
    (FDR >= !!pval) | (abs(logFC) < !!fc) ~ "NotDE"
  ))

changes <- lrt_tab %>%
  group_by(AbundanceBin) %>%
  summarize(up = sum(change == "Up"), down = sum(change == "Down"), nochange = sum(change == "NotDE"), total = n(), averExp = mean(logCPM)) %>%
  as.data.frame()

sink(paste0(odir, "/summary", ".", annot_group, suffix, ".txt"))
  print(changes)
  print(colSums(changes[, c("up", "down", "nochange")]))
sink()

# Reorganize and sort
lrt_tab <- lrt_tab %>%
  arrange(FDR, abs(logFC), PValue)

lrt_tab <- lrt_tab %>%
  left_join(coords %>% rename(pirna = name, name = name_merged))

# Split back selected/not selected seeds
lrt_tab <- lrt_tab[, c("name", "pirna", "pirna_seed", "gene_id", "AbundanceBin", "change", "logFC", "logCPM", "F", "PValue", "FDR", grep("_normCounts", colnames(lrt_tab), value = T), grep("_rawCounts", colnames(lrt_tab), value = T), grep("_coords", colnames(lrt_tab), value = T))]

lrt_tab_full <- lrt_tab %>%
  distinct() # Save everything, including coords
lrt_tab <- lrt_tab %>%
  select(!contains("_coords"), -pirna) %>%
  distinct() # Save only DE results and piRNA names

write.table(x = lrt_tab, file = paste0(odir, "/de", ".", annot_group, suffix, ".tsv"), sep = "\t", row.names = F, quote = F, col.names = T)
write.table(x = lrt_tab_full, file = paste0(odir, "/de", ".", annot_group, suffix, ".inclCoords.tsv"), sep = "\t", row.names = F, quote = F, col.names = T)