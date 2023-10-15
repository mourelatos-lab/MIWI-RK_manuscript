#!/usr/bin/env Rscript
#
# Plot ping-pong signature from bedtools window
#

library("optparse")
# library("rio")
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2"))

################################################################################

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "Input table from bedtools. Can be gz unless stdin. For stdin put \"stdin\".", metavar = "File"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    help = "Output pdf file.", metavar = "Prefix"
  ),
  make_option(
    c("-m", "--multi"),
    type = "logical", action = "store_true", default = FALSE,
    help = "Allow multimapped reads as well. Default: FALSE"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

################################################################################

read_min <- 24 # Min read length to consider
read_max <- 32 # max read length to consider
max_dist <- 30 # Maximum distance for the ping-pong calculation

ifile <- opt$input
multi <- opt$multi
ofile <- opt$output

if (ifile == "stdin") {
  tab <- read.table(file(ifile), sep = "\t", stringsAsFactors = F, header = T)
} else {
  if (tools::file_ext(ifile) == "gz") {
    tab <- read.table(gzfile(ifile), sep = "\t", stringsAsFactors = F, header = T)
  } else {
    tab <- rio::import(ifile, format = "tsv")
  }
}

if (!multi) {
  print("Keeping only uniquely mapped reads.")
  tab <- tab %>%
    filter(score == 255, score_ovl == 255)
} else {
  print("Keeping both uniquely and multi mapped reads.")
}

tab <- tab %>%
  filter(
    length >= !!read_min, length <= !!read_max,
    length_ovl >= !!read_min, length_ovl <= !!read_max
  )

# Get total number of reads and pairs
tot_pairs <- length(unique(paste(tab$name, tab$name_ovl, sep = ".")))

# Get distances
dist <- abs(tab$start - tab$start_ovl) # V2 is start at +, V10 is start at -
dist <- dist + 1 # Adjust to 0-based bed input
dist <- as.data.frame(table(dist)) %>%
  rename("Distance" = dist, "Pairs" = Freq)

dist <- dist %>%
  filter(as.numeric(Distance) <= !!max_dist)

# Get Z-score
background <- dist[dist$Distance != 10, "Pairs"]

# dist$zscore <- NA
zscore <- (dist[dist$Distance == 10, "Pairs"] - mean(background)) / sd(background) # z-score in Zamore https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3236501/; "The Ping-Pong z-score was then the difference of the score at the 52-to-52 distance of 10 nt and the mean scores of background distances, divided by the standard deviation of the scores of background distances, defined as distances of 0–9 and 11–20 nt."

p1 <- ggplot(data = dist, aes(x = Distance, y = Pairs)) +
  geom_bar(stat = "identity") +
  geom_text(data = dist %>% filter(Distance == 10), aes(x = Distance, y = Pairs, hjust = -0.25, vjust = 1, label = paste("Z-score:", round(zscore, 2)))) +
  ggtitle(paste("Total number of pairs:", tot_pairs)) +
  xlab("5'-5' Distance") +
  theme_classic()

pdf(file = ofile)
  p1
dev.off()
