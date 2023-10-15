#!/usr/bin/env Rscript
#
# Get "differential expression" of piRNA positions
# We divide the position according to the annotation - piRNA clusters, repeats, genes, intergenic
#   genes are further divided into cds, 3utr, 5utr; repeats to - LINE, SINE, other
#
################################################################################
### annotations - BED format
### IMPORTANT: These are hardcoded - please adjust to your needs/structure
## piRNA clusters
annot_clus <- "/home/joppelt/projects/pirna_mouse/data/mm10/piRNAclusterDB/pirna-clusters.bed" # unstranded!
## Repeats
annot_repeat <- "/home/joppelt/projects/pirna_mouse/data/mm10/rmsk_categ.bed" # remove rRNA, tRNA, snRNA, srpRNA, RNA
## Exon = exons including UTRs (=genes)
annot_exon <- "/home/joppelt/projects/pirna_mouse/data/mm10/ensembl_genes.exon.bed"
## CDS = exons without UTRs
annot_cds <- "/home/joppelt/projects/pirna_mouse/data/mm10/ensembl_genes.cds.bed"
## 5UTR
annot_fiveutr <- "/home/joppelt/projects/pirna_mouse/data/mm10/ensembl_genes.five_prime_utr.bed"
## 3UTR
annot_threeutr <- "/home/joppelt/projects/pirna_mouse/data/mm10/ensembl_genes.three_prime_utr.bed"
## Introns
annot_intron <- "/home/joppelt/projects/pirna_mouse/data/mm10/ensembl_genes.intron.bed"
## Intragenic =
annot_intra <- "/home/joppelt/projects/pirna_mouse/data/mm10/ensembl_genes.intergenic.bed"

################################################################################

suppressPackageStartupMessages(library("optparse"))
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

annot_granges <- function(granges, annot, annot_name, ignore.strand = TRUE, annot_field) {
  # Note: Run unique() on the resulting annotated object to remove duplicates. This is ON PURPOSE not included here
  #     in case you want to keep duplicates.
  if (length(names(mcols(granges))) > 0 && names(mcols(granges)) %in% c(annot_name, paste(annot_field, "y", sep = "."))) {
    stop(paste0("The new annotation name ", annot_name, " or ", annot_field, ".y temp version already exists in the input. Please rename or remove it and rerun."))
  }

  names_keep <- c(names(mcols(granges)), paste(annot_field, "y", sep = ".")) # Keep original columns + the one we want to add
  names(mcols(annot)) <- paste(names(mcols(annot)), "y", sep = ".") # Append ".y" to the added columns

  # Get overlap
  if (!ignore.strand) { # stranded
    granges <- plyranges::join_overlap_left_directed(granges, annot)
  } else { # unstranded
    granges <- plyranges::join_overlap_left(granges, annot)
  }

  # Keep only original columns + the newly annotated
  mcols(granges) <- mcols(granges)[names(mcols(granges)) %in% names_keep]
  names_ind <- names(mcols(granges)) == paste(annot_field, "y", sep = ".")
  names(mcols(granges))[names_ind] <- annot_name

  return(granges)
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
    c("-n", "--top"),
    type = "integer", default = 30,
    help = "Name top N genes in the output plos. Not used right now. Default: 30."
  ),
  make_option(
    c("-t", "--type"),
    type = "character", default = "pos",
    help = "Summarized counts by positions or by sequence [pos|seq]. Default: pos"
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

if (type == "pos") {
  coords <- tab_main %>%
    select(chr, start, end, name, strand) # Save coords for later
} else {
  # TODO: Parse coords from merged (|) columns
  coords <- tab_main %>%
    select(name, paste(samples, "coords", sep = "_"))
  coords <- coords[rowSums(is.na(coords)) < ncol(coords) - 1, ] # Remove sequences not present in any of the selected samples
  coords$coords <- apply(coords[, 2:ncol(coords)], 1, paste, collapse = "|")
  coords <- coords %>% select(-ends_with("_coords"))
  coords_l <- sapply(coords[, "coords"], function(x) strsplit(x = x, split = "|", fixed = T))
  coords_l <- sapply(coords_l, unique) # get only unique coords
  coords_l <- sapply(coords_l, function(x) x[x != "NA"]) # remove NA values
  coords_l <- sapply(coords_l, function(x) x[x != "."]) # remove "." values
  names(coords_l) <- coords$name
  coords <- as.data.frame(unlist(coords_l, recursive = F, use.names = T), stringsAsFactors = FALSE)
  coords$name <- rownames(coords)
  #  coords$name <- gsub('[[:digit:]]+', '', rownames(coords)) # Remove number from name
  coords <- cbind(coords, t(sapply(coords[, 1], function(x) unlist(strsplit(x, ","))))) # Get coords
  colnames(coords)[(ncol(coords) - 3):ncol(coords)] <- c("chr", "start", "end", "strand")
  coords <- coords %>%
    select(chr, start, end, name, strand)
}

tab_main <- tab_main %>%
  select(name, !!samples) %>%
  replace(is.na(.), 0)
rownames(tab_main) <- tab_main$name
tab_main$name <- NULL

# Reorder count file
tab_main <- tab_main %>%
  select(!!coldata$sample)

# Make sure coldata and counts table have the same order
if (sum((colnames(tab_main) == coldata$sample) == F) != 0) {
  stop(print("Count table and coldata table don't have the same order, please check and rerun"))
} else { # rename to a nice name
  colnames(tab_main) <- coldata$name
}

# Bring back names column
tab_main$name <- rownames(tab_main)

# Get suffix for output files
if (byfeature) {
  suffix <- "-byFeat"
} else {
  if (type == "pos") {
    suffix <- "-byPos"
  } else {
    suffix <- "-bySeq"
  }
}

# Read bed annotations
annot_clus <- read_bed(annot_clus) %>%
  coords_to_name() %>%
  bed_to_granges()
annot_repeat <- read_bed(annot_repeat) %>%
  coords_to_name() %>%
  bed_to_granges()
annot_exon <- read_bed(annot_exon) %>%
  coords_to_name() %>%
  bed_to_granges()
annot_cds <- read_bed(annot_cds) %>%
  coords_to_name() %>%
  bed_to_granges()
annot_fiveutr <- read_bed(annot_fiveutr) %>%
  coords_to_name() %>%
  bed_to_granges()
annot_threeutr <- read_bed(annot_threeutr) %>%
  coords_to_name() %>%
  bed_to_granges()
annot_intron <- read_bed(annot_intron) %>%
  coords_to_name() %>%
  bed_to_granges()
annot_intra <- read_bed(annot_intra) %>%
  coords_to_name() %>%
  bed_to_granges()

coords <- coords %>%
  bed_to_granges()

if (type == "seq") {
  coords$name <- gsub("[[:digit:]]+", "", coords$name) # Remove number from name column
}

# Annotated coordinates with annotation beds
coords <- annot_granges(coords, annot_clus, "annot.cluster", ignore.strand = TRUE, annot_field = "name") %>% unique() # same.strand = FALSE
coords <- annot_granges(coords, annot_repeat, "annot.repeat", ignore.strand = FALSE, annot_field = "name") %>% unique() # same.strand = TRUE
coords <- annot_granges(coords, annot_exon, "annot.exon", ignore.strand = FALSE, annot_field = "name") %>% unique() # same.strand = TRUE
coords <- annot_granges(coords, annot_cds, "annot.cds", ignore.strand = FALSE, annot_field = "name") %>% unique() # same.strand = TRUE
coords <- annot_granges(coords, annot_fiveutr, "annot.5utr", ignore.strand = FALSE, annot_field = "name") %>% unique() # same.strand = TRUE
coords <- annot_granges(coords, annot_threeutr, "annot.3utr", ignore.strand = FALSE, annot_field = "name") %>% unique() # same.strand = TRUE
coords <- annot_granges(coords, annot_intron, "annot.intron", ignore.strand = FALSE, annot_field = "name") %>% unique() # same.strand = TRUE
coords <- annot_granges(coords, annot_intra, "annot.intragenic", ignore.strand = FALSE, annot_field = "name") %>% unique() # same.strand = TRUE
# This should be 0 if we have all the annotations
# sum(rowSums(is.na(elementMetadata(coords)))==(ncol(elementMetadata(coords))-1))

# Extract only annotation with names
annots <- as.data.frame(elementMetadata(coords))
annots[is.na(annots)] <- "." # Replace NAs to make our lives easier

annots <- annots %>%
  distinct()

################################################################################
colkeep <- colnames(tab_main) # For easier filtering later

tab_main <- tab_main %>%
  left_join(annots)

# Run DE for each group separately
for (annot_group in c("all", colnames(annots)[!colnames(annots) %in% c("name")])) {
  if (annot_group == "all" & !byfeature) { # if this is a first run, make results for all the sequences/positions regardless on feature/group
    # if this is a first run, make results for all the sequences/positions regardless on feature/group
    print(annot_group)
    tab <- tab_main %>%
      select(!!colkeep) %>%
      distinct()
  } else if (annot_group == "all" & byfeature) {
    # We don't want all and by feature, skip
    next
  } else {
    print(annot_group)
    tab <- tab_main %>%
      filter(name %in% annots[annots[, annot_group] != ".", "name"]) %>%
      select(!!colkeep, !!annot_group) %>%
      filter(!!as.symbol(annot_group) != ".") %>% # we might have multiple lines due to multi-annotation "issue", it's safe to remove them
      distinct()

    # If we analyze by feature instead of input lines
    if (byfeature) {
      print("Summarizing input by annotation table features.")
      tab_feat <- reshape2::melt(tab, id.vars = c(annot_group, "name"))
      tab_feat <- tab_feat %>%
        group_by(variable, !!as.symbol(annot_group)) %>%
        mutate(count_feat = sum(value)) %>%
        select(variable, !!as.symbol(annot_group), count_feat) %>%
        distinct()
      tab_feat <- reshape2::dcast(tab_feat, as.formula(sprintf("%s ~ variable", annot_group)), value.var = "count_feat") %>% # We have to use "as.formula" if we want to use variable in dcast
        rename(name = !!annot_group)
      tab <- tab_feat
    } else {
      tab[, annot_group] <- NULL
      tab <- tab %>%
        distinct()
    }
  }
  # Make sure it's expressed in at least one sample
  rownames(tab) <- tab$name
  tab$name <- NULL
  tab <- tab[rowSums(tab) != 0, ]

  ### "Differential expression"
  d <- DGEList(counts = tab, group = coldata$condition) # edgeR DGE object

  design <- model.matrix(~ 0 + batch + condition, data = coldata)

  keep <- filterByExpr(d, design)
  d <- d[keep, , keep.lib.sizes = FALSE]
  d <- calcNormFactors(d) # Calculate normalization factors

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

  # # With batch effect removal
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
  lrt_tab <- lrt_tab[, c("name", "AbundanceBin", "change", "logFC", "logCPM", "F", "PValue", "FDR", grep("_normCounts", colnames(lrt_tab), value = T), grep("_rawCounts", colnames(lrt_tab), value = T))]

  write.table(x = lrt_tab, file = paste0(odir, "/de", ".", annot_group, suffix, ".tsv"), sep = "\t", row.names = F, quote = F)
}