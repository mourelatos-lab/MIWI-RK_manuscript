#!/usr/bin/env Rscript
#
# Get piRNA counts from bedtools bamtobed and compare them and their distribution
#
# Two arguments: input table in bed format; and keep only unique sequence [TRUE;FALSE] (optional; default: FALSE)
# Input can be stdin (use 'stdin')
#
################################################################################

suppressPackageStartupMessages(library("optparse"))
library("ggplot2")
suppressPackageStartupMessages(library("GenomicRanges"))
library("stringr")
suppressPackageStartupMessages(library("dplyr"))

################################################################################

fix_coords_starts <- function(bed) { # Fix coordinates if we have only start position for +/- strand regions
  bed <- bed %>%
    mutate(start = if_else(strand == "+", start, end - 1L)) %>%
    mutate(end = if_else(strand == "+", start + 1L, end))
  return(bed)
}
fix_coords_ends <- function(bed) { # Fix coordinates if we have only end position for +/- strand regions
  bed <- bed %>%
    mutate(end = if_else(strand == "+", end, start + 1L)) %>%
    mutate(start = if_else(strand == "+", end - 1L, start))
  return(bed)
}

################################################################################

option_list <- list(
  make_option(
    c("-i", "--ifile"),
    type = "character",
    help = "Input bed file from bedtools bamtobed with a header. Can be stdin (use 'stdin').", metavar = "File"
  ),
  make_option(
    c("-o", "--ofile"),
    type = "character", default = NULL,
    help = "Output file prefix if ifile is stdin. Otherwise uses the same prefix as ifile.", metavar = "File"
  ),
  make_option(
    c("-u", "--unique"),
    type = "logical", default = FALSE, action = "store_true",
    help = "Use only uniquelly mapped reads (-q 255 as in STAR). Default: all reads.", metavar = "Boolean"
  ),
  make_option(
    c("-t", "--type"),
    type = "character", default = "pos",
    help = "Summarize the counts by positions or by sequence [pos|seq]. Default: pos"
  ),
  make_option(
    c("-g", "--genome"),
    type = "logical", default = FALSE, action = "store_true",
    help = "If --type==\"seq\" does the input contain genomic sequence? Expects to see \"seq_genome\" column. If not, expects \"seq_reads\" column. Default: FALSE"
  ),
  make_option(
    c("-s", "--trimseq"),
    type = "integer", default = NULL,
    help = "If type==\"seq\" trim the sequences to maximum of trimseq. For example, for piRNAs it could be 25. Default: do not trim."
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

ifile <- opt$ifile
uniq <- opt$unique
type <- opt$type
trim_len <- opt$trimseq
genome <- opt$genome

################################################################################

if (uniq) {
  suffix1 <- "uniq"
} else {
  suffix1 <- "all"
}

if (ifile == "stdin") {
  if (is.null(opt$ofile)) {
    stop("If you use \'stdin\' you have to specify prefix/basename of the output file.")
  }
  intab <- read.table(file("stdin"), sep = "\t", header = T, stringsAsFactors = F, row.names = NULL)
} else {
  intab <- rio::import(ifile, format = "tsv")
}

if (uniq) {
  print("Getting only unique alignments...")
  before <- nrow(intab)
  intab <- intab %>%
    filter(score == 255)
  print(paste0("Kept: ", round(nrow(intab) / (before / 100), 2), "% of the input."))
} else {
  print("Using all alignments.")
}

# Get piRNA len
intab$len <- intab$end - intab$start

## Make new "name" according to mapping position
# Full
intab$full_pos <- paste(intab$chr, ",", intab$start, ",", intab$end, ",", intab$strand, sep = "")

# Only start positions
intab <- intab %>%
  mutate(
    starts = if_else(strand == "+", paste(intab$chr, ",", intab$start, ",", intab$strand, sep = ""),
      paste(intab$chr, ",", intab$end, ",", intab$strand, sep = "")
    ),
    ends = if_else(strand == "+", paste(intab$chr, ",", intab$end, ",", intab$strand, sep = ""),
      paste(intab$chr, ",", intab$start, ",", intab$strand, sep = "")
    )
  )

counts <- list()
if (type == "pos") {
  print("Counting by positions.")

  ## Get counts for full position, starts and ends
  # Make temp column in case we are merging seqeunces later
  counts$full_pos <- as.data.frame(table(intab$full_pos), stringsAsFactors = F) %>%
    rename(full_pos = Var1) %>%
    left_join(intab %>% select(-starts, -ends, -name, -score, -len) %>% distinct()) %>%
    rename(name = full_pos, score = Freq) %>%
    select(chr, start, end, name, score, strand) %>%
    distinct()
  counts$starts <- as.data.frame(table(intab$starts), stringsAsFactors = F) %>%
    rename(starts = Var1) %>%
    left_join(intab %>% select(-full_pos, -ends, -name, -score, -len) %>% distinct()) %>%
    rename(name = starts, score = Freq) %>%
    select(chr, start, end, name, score, strand) %>%
    fix_coords_starts() %>%
    distinct()
  counts$ends <- as.data.frame(table(intab$ends), stringsAsFactors = F) %>%
    rename(ends = Var1) %>%
    left_join(intab %>% select(-full_pos, -starts, -name, -score, -len) %>% distinct()) %>%
    rename(name = ends, score = Freq) %>%
    select(chr, start, end, name, score, strand) %>%
    fix_coords_ends() %>%
    distinct()

  suffix2 <- "bed"
} else if (type == "seq") {
  print("Counting by sequences")

  if (!is.null(trim_len)) {
    if (genome) {
      intab$seq <- sapply(intab$seq_genome, function(x) str_trunc(x, trim_len, "right", ellipsis = "")) # Trim sequence to trim_len
    } else {
      intab$seq <- sapply(intab$seq_reads, function(x) str_trunc(x, trim_len, "right", ellipsis = "")) # Trim sequence to trim_len
    }
    suffix0 <- paste0("trim", trim_len)
  } else {
    if (genome) {
      intab$seq <- intab$seq_genome
    } else {
      intab$seq <- intab$seq_reads
    }
    suffix0 <- "complete"
  }

  # Collapse all positions to all the unique sequences
  counts$full_pos <- intab %>%
    group_by(seq) %>%
    summarize(
      count = n(),
      full_pos = full_pos %>% unique() %>% paste(collapse = "|"),
      starts = starts %>% unique() %>% paste(collapse = "|"),
      ends = ends %>% unique() %>% paste(collapse = "|"),
      seq_reads = seq_reads %>% unique() %>% paste(collapse = "|")
    ) %>%
    rename(name = seq)

  suffix2 <- "tsv"
}

# Get the first nucleotide of the sequences
# This is very approximate because we might have collapsed sequences at the same positions
if (type == "seq") {
  nucl <- counts$full_pos
  nucl$firstnucl <- sapply(nucl$name, function(x) substr(x, 1, 1))
  nucl <- nucl %>%
    group_by(firstnucl) %>%
    summarize(count = sum(count)) %>%
    mutate(perc = round(count / (sum(count) / 100), 2))

  p0 <- ggplot(nucl, aes(x = firstnucl, y = perc, fill = firstnucl)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw()

  pdf(paste0(opt$ofile, ".pirna-firstnucl.", suffix1, ".pdf"))
    print(p0)
  dev.off()
}

if (ifile == "stdin") {
  if (type == "seq") {
    rio::export(counts$full_pos, file = paste0(opt$ofile, ".counts.", suffix0, ".", suffix1, ".", suffix2), format = "tsv")
  }
  if (type == "pos") {
    rio::export(counts$full_pos, file = paste0(opt$ofile, ".counts-full.", suffix1, ".", suffix2), format = "tsv")
    rio::export(counts$starts, file = paste0(opt$ofile, ".counts-starts.", suffix1, ".", suffix2), format = "tsv")
    rio::export(counts$ends, file = paste0(opt$ofile, ".counts-ends.", suffix1, ".", suffix2), format = "tsv")
  }
} else {
  if (type == "seq") {
    rio::export(counts$full_pos, file = gsub(pattern = ".bed", ".counts", suffix0, ".", suffix1, ".", suffix2, x = ifile), format = "tsv")
  }
  if (type == "pos") {
    rio::export(counts$full_pos, file = gsub(pattern = ".bed", ".counts-full.", suffix1, ".", suffix2, x = ifile), format = "tsv")
    rio::export(counts$starts, file = gsub(pattern = ".bed", ".counts-starts.", suffix1, ".", suffix2, x = ifile), format = "tsv")
    rio::export(counts$ends, file = gsub(pattern = ".bed", ".counts-ends.", suffix1, ".", suffix2, x = ifile), format = "tsv")
  }
}
