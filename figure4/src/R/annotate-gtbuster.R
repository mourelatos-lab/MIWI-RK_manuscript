#!/usr/bin/env Rscript
#
# Annotate file containing "flat" bed coordinates in one of the columns with gtf,
#   summarize (collapse) them by piRNAs results, and  make them more compact
# Can annotate to genes or to transcripts
# Can use GTF (Ensembl) or BED as input annotation
#
################################################################################

# library("rio")
suppressPackageStartupMessages(library("optparse"))
library("ggplot2")
suppressPackageStartupMessages(library("rtracklayer"))
library("GenomicRanges")
suppressPackageStartupMessages(library("dplyr"))
# library("plyranges")

################################################################################

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

# "Unfold" name column which has names merged over "_" for multiple positions within replicate and "|" between replicates
# We want one coorinate per row to simplify the annotation
expand_by_col <- function(tab, col, sep) {
  tab <- tab %>%
    tidyr::separate_rows(!!col, sep = as.symbol(sep), convert = FALSE) %>%
    as.data.frame()
  return(tab)
}

################################################################################

option_list <- list(
  make_option(
    c("-i", "--ifile"),
    type = "character",
    help = "Input table from GTBuster. For example: gt2.filt.tab ", metavar = "File"
  ),
  make_option(
    c("-a", "--annot"),
    type = "character",
    help = "Gene/transcript annotation in GTF or BED", metavar = "File"
  ),
  make_option(
    c("-o", "--ofile"),
    type = "character",
    help = "Output file.", metavar = "File"
  ),
  make_option(
    c("-r", "--region"),
    type = "character", default = "all",
    help = "Gene/transcript regions to use for annotation. For example: \"three_prime_utr,five_prime_utr,CDS\". Has to be split by comma if multiple. Default: use all regions."
  ),
  make_option(
    c("-d", "--database"),
    type = "character", default = "other",
    help = "[ucsc|squire|other] \"Source\" of the annotation. Used to rename columns. Default: other."
  ),
  make_option(
    c("-f", "--feature"),
    type = "character", default = "gene",
    help = "[gene|transcript] What do we want to annotate to. Default: gene ", metavar = "File"
  ),
  make_option(
    c("-u", "--keepunannot"),
    type = "logical", action = "store_true", default = FALSE,
    help = "Keep unannotated features. Default: FALSE"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

intab <- opt$ifile
otab <- opt$ofile
inannot <- opt$annot
features <- opt$region
database <- opt$database
col <- 1 # name or number of a column to use for annotation (coordinates)
sep <- "," # Separator of the column for coordinates in the input table; default: ","
center <- TRUE # If the to-be-annotated BED is a region, do you want to center the region and annotated that instead of the whole region
annot_to <- opt$feature

################################################################################

if (any(grepl(",", features))) { # comma separate features you want to analyze/subset from the annotation
  features <- unlist(strsplit(features, ","))
}

if (annot_to == "transcript") {
  annot_col <- "transcript_id"
} else {
  annot_col <- "gene_id"
}

# Read annotation gtf/bed
# if gtf
informat <- gsub(".*\\.", "", inannot)

if (informat %in% c("gtf", "bed")) {
  annot <- rtracklayer::import(inannot, format = informat)
} else {
  stop(paste0("Cannot recognize annotation input format: ", informat, ". Please check if it ends as .gtf/.bed and rerun. Or test the new format and add it to the list of allowed/tested formats."))
}

if (informat == "gtf") {
  annot$type <- as.character(annot$type) # Change from factor to character
} else if (informat == "bed") {
  if (database != "squire") {
    annot$type <- annot$name
    annot$name <- paste0(seqnames(annot), ":", start(annot), "-", end(annot), ",", strand(annot), ",", annot$name)
    mcols(annot)[, annot_col] <- annot$name
  } else {
    annot$type <- unlist(lapply(annot$name, function(x) unlist(strsplit(x, "|", fixed = T))[4]))
    mcols(annot)[, annot_col] <- annot$name
  }
}

# Subset annotation features if specified
if (length(features) != 0 && features != "all") {
  print(paste("Subsetting the annotation only for:", paste(features, collapse = ",")))
  annot <- annot[annot$type %in% features] # Keep only "informative" parts of annot
} else {
  print("Using all annotation features.")
}
# keep only relevant columns from the annotation
mcols(annot) <- mcols(annot)[, c(annot_col, "type")]
annot <- unique(annot)

# Read input table from GTBuster (filtered)
tab <- rio::import(intab, format = "tsv")

# Expand multi-coords into rows and remove expression values
tab <- expand_by_col(tab = tab, col = "V1", sep = "_|\\|")
tab$V1 <- sub(",[^,]+$", "", tab$V1)
tab <- unique(tab)

tab <- cbind(tab, t(sapply(tab[, col], function(x) unlist(strsplit(x, sep))))) # Get coords
colnames(tab)[(ncol(tab) - 3):ncol(tab)] <- c("chr", "start", "end", "strand")

tab <- tab %>%
  dplyr::rename(degra_coord = V1, pirna_pos = V2, pirna_seed = V3, pirna = V4)

tab_granges <- tab %>%
  dplyr::select(chr, start, end, strand, degra_coord, pirna_seed, pirna, pirna_pos) %>%
  bed_to_granges()

# IMPORTANT: This is specific for our degradome where we have region of degradome cut and then extra 50 bp on each side
# Extract the a single point in the middle of the range
if (center) {
  print("Centering of the region is set to TRUE. Using only midpoint position for the annotation.")
  tab_granges$start.bckp <- start(tab_granges)
  tab_granges$end.bckp <- end(tab_granges)
  tab_granges$width.bckp <- width(tab_granges)
  start(tab_granges) <- start(tab_granges) + floor(tab_granges$width.bckp / 2)
  end(tab_granges) <- end(tab_granges) - floor(tab_granges$width.bckp / 2)
}

tab_granges <- annot_granges(granges = tab_granges, annot = annot, annot_name = annot_col, ignore.strand = FALSE, annot_field = annot_col)
tab_granges <- annot_granges(granges = tab_granges, annot = annot, annot_name = "type", ignore.strand = FALSE, annot_field = "type")
if (!opt$keepunannot) {
  tab_granges <- tab_granges[!is.na(elementMetadata(tab_granges)[, annot_col]), ] # Remove unannotated
}
tab_granges <- unique(tab_granges)

# Get summary of number of unique piRNAs and piRNA seed sequences per gene
tab_granges <- tab_granges %>%
  as.data.frame() %>%
  group_by(!!sym(annot_col)) %>%
  mutate(pirna_mols = length(unique(pirna)), pirna_seeds = length(unique(pirna_seed)))

# Plot of number of targeting piRNAs per gene
p1 <- tab_granges %>%
  select(!!sym(annot_col), pirna_mols) %>%
  distinct() %>%
  ggplot(aes(pirna_mols)) +
  geom_histogram(binwidth = 1) +
  coord_cartesian(xlim = c(1, 100)) +
  theme_bw() +
  xlab("Number of piRNA molecules") +
  ylab("Genes") +
  ggtitle("Number of targeting piRNAs per gene") +
  theme(plot.title = element_text(hjust = 0.5))

# Plot of number of targeting piRNA seeds per gene
p2 <- tab_granges %>%
  select(!!sym(annot_col), pirna_seeds) %>%
  distinct() %>%
  ggplot(aes(pirna_seeds)) +
  geom_histogram(binwidth = 1) +
  coord_cartesian(xlim = c(1, 30)) +
  theme_bw() +
  xlab("Number of piRNA seeds") +
  ylab("Genes") +
  ggtitle("Number of targeting seeds per gene") +
  theme(plot.title = element_text(hjust = 0.5))

# Get "most hit" regions per gene and in general
region_per_genes <- tab_granges %>%
  group_by(!!sym(annot_col)) %>%
  count(type)

print("Total piRNA target regions")
tab_granges %>%
  group_by(!!sym(annot_col)) %>%
  count(type) %>%
  group_by(type) %>%
  tally(n)

# Get "nubmer of hits" per pirna
pirna_targets <- tab_granges %>%
  group_by(pirna) %>%
  summarize(gene_count = length(unique(!!sym(annot_col))))
# table(pirna_targets$gene_count)

pirna_targets_seed <- tab_granges %>%
  group_by(pirna_seed) %>%
  summarize(gene_count = length(unique(!!sym(annot_col))))

# Plot number of number of genes targeted by individual piRNAs
p3 <- ggplot(data = pirna_targets, aes(gene_count)) +
  geom_histogram(binwidth = 1) +
  coord_cartesian(xlim = c(1, 30)) +
  theme_bw() +
  xlab("Number of genes") +
  ylab("piRNA molecules") +
  ggtitle("Number of targeted genes per piRNA") +
  theme(plot.title = element_text(hjust = 0.5))

# Plot number of number of genes targeted by individual seeds
p4 <- ggplot(data = pirna_targets_seed, aes(gene_count)) +
  geom_histogram(binwidth = 1) +
  coord_cartesian(xlim = c(1, 50)) +
  theme_bw() +
  xlab("Number of genes") +
  ylab("piRNA seeds") +
  ggtitle("Number of targeted genes per seed") +
  theme(plot.title = element_text(hjust = 0.5))

# Bring back original coords for better compatibility later on
tab_granges$start <- tab_granges$start.bckp
tab_granges$end <- tab_granges$end.bckp
tab_granges$width <- tab_granges$width.bckp
tab_granges <- tab_granges %>%
  select(-start.bckp, -end.bckp, -width.bckp)
tab_granges <- unique(tab_granges)

rio::export(x = tab_granges, file = otab, format = "tsv")

pdf(gsub(pattern = ".tab", replacement = ".pdf", x = otab))
  print(p1)
  print(p2)
  print(p3)
  print(p4)
dev.off()

### GTBuster output https://github.com/weng-lab/GTBuster
# Column 1: description of the cleavage site from degradome-seq: chr|start position|end position|strand|normalized degradome read counts under condition 1 (WT)|normalized degradome read counts under condition 2 (mut). The extra strand in the end is not useful, just a byproduct of some function for sanity check.
# Column 2: the position of the cleavage (0-based). The reference is 50nt upstream + degradome cleavage site (1nt) + 50nt downstream. 58 means the cut happens with the exactly 10nt overlap of piRNA and target sites. I believe only rows with 58 in this column are kept in these files.
# Column 3: 2–7 nt of the piRNA, i.e., the seed portion (or, at least, this is what I use for initial identification of potential target pairs)
# Column 4: piRNA sequence
# Column 5: normalized piRNA abundance (RPM)
# Column 6: the remaining portion of the piRNA (column 3 subtracted from column 4), i.e., so-called non-seed portion of piRNA (assuming 2–7 nt as the seed region)
# Column 7: the part of the target that corresponds to the non-seed part of piRNA. This is reverse-complemented to make it easier for a human to read.
# Column 8: the part of the target corresponding to the first nt of piRNA. Pairing of this position is not required, just like in Wei's paper
# Columns 9&10: please ignore these, as these are just short versions of columns 7&8
# Column 11–14: # matches, # GU wobble pairs, matching positions, GU wobble positions
