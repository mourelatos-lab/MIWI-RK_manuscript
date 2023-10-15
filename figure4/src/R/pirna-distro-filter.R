#!/usr/bin/env Rscript
#
# Basic script to compare/prefilter piRNA libraries
#
# If the compared libraries are biological replicates (filtcommon=TRUE), it keeps only piRNAs present
#   in both biological replicates and piRNAs with min.count
#
# If the compare libraries are not biological replicates (filtcommon=FALSE, default) or if the library is without
#   replicates it keep only piRNAs with min.count
#
################################################################################

suppressPackageStartupMessages(library("optparse"))
# library("stringr")
suppressPackageStartupMessages(library("dplyr"))
library("ggplot2")

################################################################################

read_bed <- function(input) {
  tab <- read.table(input, header = F, sep = "\t", stringsAsFactors = F)
  colnames(tab)[1:6] <- c("chr", "start", "end", "name", "score", "strand")
  return(tab)
}
get_samples <- function(name_vect, sepa) {
  samp_names.tmp <- t(as.data.frame(strsplit(name_vect, sepa, fixed = T)))
  ind <- lapply(apply(samp_names.tmp, 2, unique), length) != 1
  samp_names.tmp <- as.data.frame(samp_names.tmp[, ind])
  samp_names <- as.character(interaction(samp_names.tmp, sep = sepa))
  return(samp_names)
}
split_str <- function(vect, delim, pos) {
  unlist(lapply(vect, function(x) strsplit(x, split = delim, fixed = T)[[1]][pos]))
}

subtract_regions <- function(inbed, subtractbed) {
  inbedplus <- inbed %>%
    filter(strand == "+")
  inbedminus <- inbed %>%
    filter(strand == "-")
  subtractplus <- subtractbed %>%
    filter(strand == "+")
  subtractminus <- subtractbed %>%
    filter(strand == "-")
  plus <- bedr.subtract.region(x = inbedplus, y = subtractplus, check.chr = FALSE, check.merge = FALSE)
  minus <- bedr.subtract.region(x = inbedminus, y = subtractminus, check.chr = FALSE, check.merge = FALSE)
  outbed <- bedr.sort.region(rbind(plus, minus), check.chr = F, check.merge = F)
  return(outbed)
}

################################################################################

option_list <- list(
  make_option(
    c("-i", "--ifile1"),
    type = "character",
    help = "Input bed/sequence tsv - file 1", metavar = "File"
  ),
  make_option(
    c("-f", "--ifile2"),
    type = "character", default = NULL,
    help = "Input bed/sequence tsv - file 2", metavar = "File"
  ),
  make_option(
    c("-o", "--odir"),
    type = "character",
    help = "Output directory where all the plots will be saved.", metavar = "Directory"
  ),
  make_option(
    c("-m", "--mincount"),
    type = "integer", default = 10,
    help = "Minimum piRNA read count (raw) in at least one sample. All samples must have expression of at least 1!. Default: 10"
  ),
  make_option(
    c("-r", "--filtcommon"),
    type = "logical", default = "FALSE",
    help = "Keep only common piRNAs. Default: FALSE", metavar = "Boolean"
  ),
  make_option(
    c("-q", "--quant"),
    type = "double", default = 0.1,
    help = "Remove the percentage of least expressed genes [0-1]. Only used for visualization!. Default: 0.1"
  ),
  make_option(
    c("-t", "--type"),
    type = "character", default = "pos",
    help = "Summarize the counts by positions or by sequence [pos|seq]. Default: pos"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

################################################################################

ifile1 <- opt$ifile1
isifile2 <- !is.null(opt$ifile2) # Save if we have/don't have ifile2 for easier conditions
if (isifile2) {
  ifile2 <- opt$ifile2
}
odir <- opt$odir
min_count <- opt$mincount
filt_common <- opt$filtcommon # Are those biol. replicates? If TRUE we'll keep only piRNAs present in both samples
remove_quant <- opt$quant # Remove least 10% expressed genes (average of replicates)
type <- opt$type

################################################################################

if (filt_common) {
  if (!isifile2) {
    stop(print("If you have filt_common=TRUE you need to supply ifile2. Exiting."))
  }
}

if (isifile2) {
  samps <- get_samples(c(ifile1, ifile2), sepa = "/")
  samp1 <- samps[1]
  samp2 <- samps[2]
} else {
  samp1 <- t(as.data.frame(strsplit(ifile1, "/", fixed = T)))[2]
}

intab1 <- rio::import(ifile1, format = "tsv", header = T)
if (isifile2) {
  intab2 <- rio::import(ifile2, format = "tsv", header = T)
  intab <- rbind(intab1, intab2)
} else {
  intab <- intab1
}

if (type == "pos") {
  intab <- intab %>%
    rename(count = score)
}

# Save coords for later
if (type == "seq") {
  # We will dcast by sequence and keep mapped positions in all the samples
  coords <- reshape2::dcast(intab, name ~ sample, value.var = "full_pos") %>%
    replace(is.na(.), ".") %>%
    distinct()
  colnames(coords)[2:ncol(coords)] <- paste0(colnames(coords)[2:ncol(coords)], "_coords")
} else {
  coords <- intab %>%
    select(chr, start, end, name, strand) %>%
    distinct()
}

# Get all seqs with at least min_count read in at least one of the samples and more than 0 in all samples
keepme <- intab %>%
  filter(count > 0) %>%
  group_by(name) %>%
  filter(sum(count >= !!min_count) >= 1) %>% # , keep2 = sum(count!=0)
  ungroup() %>%
  pull(name) %>%
  unique()

# Get RPM and remove super-low count molecules
intab <- intab %>%
  group_by(sample) %>%
  mutate(rpm = (count) / (sum(count) / 1000000)) %>%
  ungroup() %>%
  filter(name %in% keepme)

counts_raw <- intab # Make a backup of raw counts for the future uses

intab <- reshape2::dcast(intab, formula = name ~ sample, value.var = "rpm")
intab[is.na(intab)] <- 0

# Keep only piRNAs in both replicates
if (isifile2) {
  if (filt_common) {
    print("Keeping only piRNAs common in all samples.")
    cts <- nrow(intab)
    intab <- intab[rowSums(intab[, samps] != 0) == length(samps), ] # Must be present in both samples
    print(paste0("Kept ", round(nrow(intab) / (cts / 100), 2), "% of the original number of sequences."))
  }
}

keep <- intab$name # Save piRNA names to keep (we'll do min.expression filtering in DE script)

dir.create(odir, recursive = T)

# Get the first nucleotide of the sequences
# This is very approximate because we might have collapsed sequences at the same positions
if (type == "seq") {
  counts_raw$firstnucl <- sapply(counts_raw$name, function(x) substr(x, 1, 1))
  nucl <- counts_raw %>%
    group_by(sample, firstnucl) %>%
    summarize(count = sum(count)) %>%
    group_by(sample) %>%
    mutate(perc = round(count / (sum(count) / 100), 2))

  p0 <- ggplot(nucl, aes(x = firstnucl, y = perc, fill = sample)) +
    geom_bar(stat = "identity", position = "dodge")

  if (isifile2) {
    pdf(paste0(odir, "/", samp1, "-", samp2, ".pirna-firstnucl.pdf"))
  } else {
    pdf(paste0(odir, "/", "pirna-firstnucl.pdf"))
  }
  print(p0)
  dev.off()
}

if (isifile2) { # Compare only if we have anything to compare
  # Get mean expression
  intab$meanExp <- rowSums(intab[, samps]) / length(samps)
  # Remove least x% expressed genes
  lim <- quantile(intab$meanExp, probs = remove_quant)[1]
  intab <- intab %>%
    filter(meanExp > !!lim)

  intab$logfc <- log2((intab[, samp2] + 0.01) / (intab[, samp1] + 0.01))

  intab$dir <- "nochange"
  intab$dir[intab$logfc >= 1] <- "up"
  intab$dir[intab$logfc <= (-1)] <- "down"

  # up/down/nochange by quantiles
  quants <- quantile(log2(intab$meanExp + 0.01), probs = seq(0, 1, 0.05)) # just print the quartiles
  if (any(duplicated(quants))) { # If we have any duplicated quants just add a tiny number
    print("There are duplicated quantile limits, adding a tiny number.")
    quants[duplicated(quants)] <- quants[duplicated(quants)] + runif(sum(duplicated(quants)), 1e-9, 1e-8)
  }
  intab <- within(intab, AbundanceBin <- as.integer(cut(log2(intab$meanExp + 0.01), quants, include.lowest = TRUE), include.lowest = TRUE))
  changes <- intab %>%
    group_by(AbundanceBin) %>%
    summarize(up = sum(dir %in% "up"), down = sum(dir %in% "down"), nochange = sum(dir %in% "nochange"), total = n(), averExp = mean(log2(meanExp + 0.01))) %>%
    as.data.frame()
  print(changes)
  print(colSums(changes[, c("up", "down", "nochange")]))

  if (type == "seq") {
    intab$firstnucl <- sapply(intab$name, function(x) substr(x, 1, 1))
    intab$firstnucl <- factor(intab$firstnucl, levels = unique(intab$firstnucl))

    p1 <- ggplot(intab, aes_string(x = samp1, y = samp2, color = "dir", shape = "firstnucl"))
  } else {
    p1 <- ggplot(intab, aes_string(x = samp1, y = samp2, color = "dir"))
  }
  p1 <- p1 +
    geom_point(size = 1.5, alpha = 0.5) +
    geom_point(alpha = 0.5) +
    theme_bw() +
    geom_abline(slope = 1, intercept = 0) +
    scale_color_manual(
      values = c("nochange" = "black", "up" = "red", "down" = "blue"),
      labels = c(
        paste0("Down", " (", sum(intab$dir == "down"), ")"),
        paste0("NotChanged", " (", sum(intab$dir == "nochange"), ")"),
        paste0("Up", " (", sum(intab$dir == "up"), ")")
      )
    ) +
    theme(legend.title = element_blank(), legend.position = "bottom")

  png(paste0(odir, "/", samp1, "-", samp2, ".pirna-compare.", type, ".png"), width = 960, height = 960)
    print(p1)
  dev.off()
}

# Export filtered raw counts
counts_raw <- reshape2::dcast(counts_raw, formula = name ~ sample, value.var = "count") # Get raw counts
counts_raw[is.na(counts_raw)] <- 0

counts_raw <- counts_raw %>%
  filter(name %in% keep) %>%
  left_join(coords, by = "name")

if (isifile2) {
  rio::export(x = counts_raw, file = paste0(odir, "/", samp1, "-", samp2, ".pirna-counts.filt.", type, ".tsv"), format = "tsv")
} else {
  rio::export(x = counts_raw, file = paste0(odir, "/", "pirna-counts.filt.", type, ".tsv"), format = "tsv")
}
