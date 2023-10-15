#!/usr/bin/env Rscript
#
# Converts genomic GTBuster tab (chromosome) to transcriptomic bed (transcripts)
#
# Works only for Ensembl GTF!
# https://bioconductor.org/packages/devel/bioc/vignettes/ensembldb/inst/doc/coordinate-mapping.html#2_Mapping_genomic_coordinates_to_transcript-relative_coordinates
#
#
################################################################################

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("ensembldb"))

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

################################################################################
option_list <- list(
  make_option(
    c("-i", "--ifile"),
    type = "character", default = "stdin",
    help = "Input table from GTBuster (genomic coords) with header. For example: gt2.filt.tab. For stdin use \"stdin\". Default: stdin.", metavar = "File"
  ),
  make_option(
    c("-g", "--gtf"),
    type = "character",
    help = "Ensembl named gtf to use for conversion. For example: Mus_musculus.GRCm38.99.gtf", metavar = "File"
  ),
  make_option(
    c("-o", "--ofile"),
    type = "character", default = NULL,
    help = "Output BED file. Default: same as input but with .bed", metavar = "File"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

if (opt$ifile == "stdin") {
  ifile <- read.table(file = file("stdin"), sep = "\t", stringsAsFactors = F, header = T, quote = "")
} else {
  ifile <- rio::import(opt$ifile, format = "tsv", header = T)
}

################################################################################

# Make it look like bed
ifile$start <- ifile$start - 1
ifile$name <- "."
ifile$score <- "."

ifile <- ifile %>%
  dplyr::rename("chr" = seqnames) %>%
  bed_to_granges()

center <- TRUE
if (center) {
  ifile$start.bckp <- start(ifile)
  ifile$end.bckp <- end(ifile)
  ifile$width.bckp <- width(ifile)
  start(ifile) <- start(ifile) + floor(ifile$width.bckp / 2)
  end(ifile) <- end(ifile) - floor(ifile$width.bckp / 2)
}

if (!file.exists(paste(opt$gtf, "sqlite", sep = "."))) {
  gtfTxDb <- ensDbFromGtf(opt$gtf, outfile = paste(opt$gtf, "sqlite", sep = "."))
} else {
  gtfTxDb <- paste(opt$gtf, "sqlite", sep = ".")
}
edb <- EnsDb(gtfTxDb)
gnm_tx <- genomeToTranscript(x = ifile, db = edb)
gnm_tx <- unlist(gnm_tx)

# TODO: Export IRanges to bed
df <- data.frame(
  chr = names(gnm_tx),
  start = start(gnm_tx) - 1,
  end = end(gnm_tx),
  name = c(rep(".", length(gnm_tx))),
  score = c(rep(".", length(gnm_tx))),
  strand = mcols(gnm_tx)$seq_strand,
  degra_coord_tmp = paste(mcols(gnm_tx)$seq_name, mcols(gnm_tx)$seq_start, mcols(gnm_tx)$seq_end, mcols(gnm_tx)$seq_strand, sep = ","),
  stringsAsFactors = F
)

# Add additional columns from the original input if we want to
mcols(ifile)$degra_coord_tmp <- paste(seqnames(ifile), ranges(ifile), ranges(ifile), strand(ifile), sep = ",")

df <- df %>%
  left_join(mcols(ifile) %>% as.data.frame() %>% dplyr::select(degra_coord, degra_coord_tmp, pirna_seed, pirna, pirna_mols, pirna_seeds), by = "degra_coord_tmp") %>%
  dplyr::select(-degra_coord_tmp) %>%
  distinct()

if (opt$ifile != "stdin" & is.null(opt$ofile)) {
  ofile <- paste(sub(".[^.]+$", "", opt$ifile), "bed", sep = ".")
  write.table(df, file = ofile, quote = F, sep = "\t", row.names = F, col.names = T)
} else {
  write.table(x = df, quote = F, sep = "\t", row.names = F, col.names = T)
}
