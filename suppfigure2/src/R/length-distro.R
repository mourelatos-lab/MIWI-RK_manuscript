#!/usr/bin/env Rscript
#
# Read length distribution and statistics from bedtools intersect table
#
# If you are using stdin use 'stdin'
#
################################################################################

library("optparse")
library("ggsci") # scientific pallettes
suppressMessages(library("ggrepel"))
# suppressMessages(library("rio"))
# suppressMessages(library("stringr"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("reshape2"))
# bit64

################################################################################

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "Input table from bedtools intersect with added Sample column. For stdin put \"stdin\"", metavar = "File"
  ),
  make_option(
    c("-o", "--ofile"),
    type = "character", default = "length-distro.pdf",
    help = "Output pdf file.", metavar = "Prefix"
  ),
  make_option(
    c("-f", "--from"),
    type = "integer", default = 24,
    help = "Highlight reads from xx nt. Default: 24", metavar = "Integer"
  ),
  make_option(
    c("-t", "--to"),
    type = "integer", default = 32,
    help = "Highlight reads to xx nt. Default: 32", metavar = "Integer"
  ),
  make_option(
    c("-b", "--by"),
    type = "integer", default = 1,
    help = "Size of x-axis ticks ", metavar = "Integer"
  ),
  make_option(
    c("-s", "--unstranded"),
    type = "logical", default = FALSE, action = "store_true",
    help = "Don't look at overlap strands (=make simple sense/antisense by the mapping, not by overlap). Default: Look at strands", metavar = "Boolean"
  ),
  make_option(
    c("-u", "--unique"),
    type = "logical", default = FALSE, action = "store_true",
    help = "Use only uniquely mapped reads. Default: Use all reads.", metavar = "Boolean"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

input_file <- opt$input
ofile <- opt$ofile
from <- opt$from
to <- opt$to
from_main <- 19 # plot only seq from this size
to_main <- 38 # plot only seq to this size
by <- opt$by
uniq_map <- as.logical(opt$unique)

################################################################################

if (input_file == "stdin") {
  tab <- read.table(file(input_file), header = T, sep = "\t", stringsAsFactors = F)
} else {
  tab <- read.table(input_file, header = T, sep = "\t", stringsAsFactors = F)
}

if (uniq_map) {
  print("Keeping only uniquely mapped reads (score 255)")
  tab <- tab %>%
    filter(score == 255)
} else {
  print("--unique not specified, using all mapped reads")
}

# Check duplicated reads and keep only one
if (sum(duplicated(tab$name)) > 0) {
  print("There are some duplicated read names (=reads), keep only one.")
  tab <- tab[!duplicated(tab$name), ]
}

if (!any(names(tab) == "Sample")) {
  tab$Sample <- basename(input_file)
}

# Get length tables for sense and antisense mapping separately
if (!as.logical(opt$unstranded) & any(names(tab) == "strand_ovl")) { # If we have results from overlaps
  tab <- tab %>%
    mutate(length = end - start, Strand = ifelse(strand %in% c("+", "-") & strand_ovl %in% c("+", "-"), ifelse(strand == strand_ovl, "sense", "antisense"), "not_determined")) %>%
    group_by(Sample, length, Strand) %>%
    summarize(reads = n())
} else {
  print("Making simple +/- strand mapping")
  tab <- tab %>%
    mutate(length = end - start, Strand = ifelse(strand %in% c("+", "-"), ifelse(strand == "+", "sense", "antisense"), "not_determined")) %>%
    group_by(Sample, length, Strand) %>%
    summarize(reads = n())
}

# Calculate percentages
tab <- tab %>%
  group_by(Sample) %>%
  mutate(perc = reads / (sum(reads) / 100))

y_max <- tab %>%
  group_by(Sample, length) %>%
  summarize(reads = sum(reads), perc = sum(perc)) %>%
  pull(perc) %>%
  max()
x_lab <- seq(from_main, to_main, by = 1)

tab$Sample <- factor(tab$Sample, levels = sort(unique(tab$Sample)))

mypal <- pal_npg("nrc", alpha = 1)(nlevels(tab$Sample))
cols <- mypal
names(cols) <- levels(tab$Sample)

### Make percentage plots (both sense and antisense together)
distr_plot <- tab %>%
  group_by(Sample, length) %>%
  summarize(reads = sum(reads), perc = sum(perc)) %>%
  ggplot(aes(x = length, y = perc)) +
  geom_line(aes(color = Sample), na.rm = TRUE) +
  geom_point(aes(color = Sample, shape = Sample), alpha = .5, size = 1.5) +
  scale_shape_manual(values = 1:nlevels(tab$Sample)) +
  scale_x_continuous(limits = c(from_main, to_main), breaks = seq(from_main, to_main, by = 1), labels = x_lab) +
  scale_y_continuous(breaks = seq(0, y_max + 2, by = 2)) +
  scale_color_manual(values = cols) +
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  theme_classic() +
  theme(legend.position = "bottom", legend.text = element_text(size = 10)) +
  ggtitle(paste("Read length distribution (capped at", to_main, "bp)")) +
  xlab("Length (nt)") +
  ylab("Percentage")

if (!is.null(from) & !is.null(to) & !is.null(by)) {
  distr_plot <- distr_plot + geom_vline(
    xintercept = seq(from, to, by = by), linetype = "dotted",
    color = "grey", size = 0.5
  )
}

### Make raw counts plots for plus and minus
tab <- tab %>%
  mutate(reads = ifelse(Strand == "not_determined", reads, ifelse(Strand == "sense", reads, -reads)))

# Get the "magnitude" of the number for y-axis brakes and y-axis breaks
log10_ceiling <- function(x) {
  10^(floor(log10(x)))
}

# Make all plots and plot them at once at the end https://stackoverflow.com/a/31994539/1393566
plot_data <- function(tab, samp) {
  tab_samp <- tab %>%
    filter(Sample == !!samp)

  # Add stats to the plot
  reads_summary <- tab_samp %>%
    group_by(Strand) %>%
    dplyr::summarize("Reads" = sum(abs(reads))) %>%
    ungroup() %>%
    mutate("Percentage" = Reads / (sum(Reads) / 100))
  reads_summary <- reads_summary[order(reads_summary$Strand), ]

  annotations <- data.frame(
    X = c(Inf),
    Y = c(Inf),
    text = paste0(reads_summary$Strand, ": ", reads_summary$Reads, " reads, ", round(reads_summary$Percentage, 2), "%", collapse = "\n"),
    x_adjust = c(1),
    y_adjust = c(1)
  )

  # Make y-axis and x-axis limits
  y_scale_break <- log10_ceiling(max(abs(tab_samp$reads)))

  y_lims <- c(min(tab_samp$reads) * 1.05, max(tab_samp$reads) * 1.05)
  y_min <- floor(min(y_lims) / (y_scale_break / 10)) * (y_scale_break / 10)
  y_max <- ceiling(max(y_lims) / (y_scale_break / 10)) * (y_scale_break / 10)

  if (y_min < 0) {
    y_scale <- seq(0, abs(y_min), by = y_scale_break)
    y_scale <- sort(-y_scale)
    y_scale <- c(y_scale, seq(0, y_max, by = y_scale_break))
  } else {
    y_scale <- seq(0, y_max, by = y_scale_break)
  }

  # Get colors
  tab_samp$Strand <- factor(tab_samp$Strand, levels = c("sense", "antisense", "not_determined")) # sort(unique(tab_samp$Strand)))

  mypal <- pal_npg("nrc", alpha = 1)(4)
  cols <- mypal[c(1, 3, 4)]
  names(cols) <- c("antisense", "not_determined", "sense") # levels(tab_samp$Strand)
  cols <- cols[names(cols) %in% unique(tab_samp$Strand)]

  length_plot <- ggplot(tab_samp, aes(length)) +
    geom_bar(
      data = subset(tab_samp, Strand == "sense"),
      aes(y = reads, fill = Strand), stat = "identity", position = "dodge"
    ) +
    geom_bar(
      data = subset(tab_samp, Strand == "antisense"),
      aes(y = reads, fill = Strand), stat = "identity", position = "dodge"
    ) +
    geom_bar(
      data = subset(tab_samp, Strand == "not_determined"),
      aes(y = reads, fill = Strand), stat = "identity", position = "dodge"
    ) +
    geom_hline(yintercept = 0, colour = "grey90") +
    scale_x_continuous(limits = c(from_main, to_main), breaks = seq(from_main, to_main, by = 1), labels = x_lab) +
    scale_y_continuous(breaks = y_scale) +
    scale_fill_manual(values = cols) +
    geom_text(data = annotations, aes(
      x = X, y = Y, hjust = x_adjust, vjust = y_adjust, label = text
    )) +
    theme_classic() +
    theme(legend.position = "bottom", legend.text = element_text(size = 10)) +
    xlab("Length (nt)") +
    ylab("Reads") +
    ggtitle(samp)

  print(length_plot)
}

myplots <- lapply(levels(tab$Sample), function(x) plot_data(tab = tab, samp = x))

# Plots
pdf(ofile, width = 10)
  print(
    distr_plot
  )
  for (p1 in 1:length(myplots)) {
    print(myplots[p1])
  }
dev.off()

write.table(x = tab, file = gsub(pattern = ".pdf", ".tsv", x = ofile), quote = F, sep = "\t", col.names = T, row.names = F)