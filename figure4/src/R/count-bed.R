#!/usr/bin/env Rscript
#
# Very rough and dirty script to get unique alignments (255 in score - STAR alignment) from
#   BED file, count them and collapse them and export them
# Assumes stdin; outputs stdout
#
################################################################################

suppressPackageStartupMessages(library("dtplyr")) # testing dtplyr backend
suppressPackageStartupMessages(library("dplyr"))

################################################################################

intab <- read.table(file = file("stdin"), sep = "\t", header = F, stringsAsFactors = F)
cols <- c("chr", "start", "end", "name", "score", "strand")
colnames(intab) <- cols

intab$name <- paste(intab$chr, intab$start, intab$end, intab$strand, sep = ",")

intab <- lazy_dt(intab)

intab <- intab %>%
  filter(score == "255") %>% # Get only unique mappings
  group_by(name) %>%
  summarize(
    chr = unique(chr),
    start = unique(start),
    end = unique(end),
    score = n(),
    strand = unique(strand)
  ) %>%
  ungroup() %>%
  mutate(rpm = score / (sum(score) / 1000000)) %>%
  as.data.frame()

intab$name <- paste(intab$name, round(intab$rpm, 3), sep = ",")

intab <- intab[, c(cols, "rpm")]

write.table(
  x = intab,
  file = "",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE, col.names = FALSE
)