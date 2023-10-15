#!/usr/bin/env Rscript
#
# Filter GTBuster results from randomly generated transcript positions WITHOUT filtering for cut position (58 in the normal script filter-gtbuster-nopos.R)
#
################################################################################

suppressPackageStartupMessages(library("dplyr"))

filter_gtbuster <- function(tab, pirna_pos = NULL) {
  # Keep only lines where piRNA binds from specific location

  if (!is.null(pirna_pos)) { # If we filter by position of piRNA relative to the analyzed fragment
    tab <- tab %>%
      filter(V2 == !!pirna_pos)
  }

  tab <- tab %>%
    filter(V11 >= !!min_match_full) %>%
    filter((V11 + V12) >= (!!min_match_full + !!match_wobble))
  tab$V1 <- sub(",[^,]+$", "", tab$V1) # Remove expression values
  return(tab)
}

args <- commandArgs(trailingOnly = TRUE)

repl1 <- args[1]
if (length(args) >= 3) {
  otab <- args[3]
} else {
  otab <- args[2]
}

# Check if we set secondary analysis otherwise use default
if (length(args) == 4) {
  secondary <- as.logical(args[4])
} else {
  secondary <- FALSE
}

################################################################################
# Combination of 14 and 0 gives 14 which allows 4 mismatches (~80% full match) in the 18 "searchable" (25 - 1 - 6) string
min_match_full <- 14 # Minimum of full matches after seed (has to have full match of 2-7) and before 25 nt (GTBuster uses only first 25 bp); total of "searchable" is 18 (25 - 1 - 6)
match_wobble <- 0 # How many wobbles we allow in total sum of matches

if (substr(repl1, nchar(repl1) - 2, nchar(repl1)) == ".gz") {
  repl1 <- read.table(file = gzfile(repl1), sep = "\t", header = F, stringsAsFactor = F)
} else {
  repl1 <- rio::import(repl1, format = "tsv")
}
print("Number of unfiltered positions in replicate/sample 1:")
length(unique(repl1$V1))

if (!secondary) {
  repl1 <- filter_gtbuster(repl1)
  print("Number of filtered positions in replicate/sample 1:")
  print(length(unique(repl1$V1)))
}

if (length(args) >= 3) {
  repl2 <- args[2]
  if (substr(repl2, nchar(repl2) - 2, nchar(repl2)) == ".gz") {
    repl2 <- read.table(file = gzfile(repl2), sep = "\t", header = F, stringsAsFactor = F)
  } else {
    repl2 <- rio::import(repl2, format = "tsv")
  }
  print("Number of unfiltered positions in replicate/sample 2:")
  print(length(unique(repl2$V1)))

  if (!secondary) {
    repl2 <- filter_gtbuster(repl2)
    print("Number of filtered positions in replicate/sample 1:")
    print(length(unique(repl2$V1)))
  }
} else {
  repl2 <- NULL
}

# Subset cuts in both replicates only
if (!is.null(repl2)) {
  print("Cuts in replicate 1, replicate 2 and common cuts for both replicates/samples:")

  repl1_cuts <- unique(repl1$V1)
  repl2_cuts <- unique(repl2$V1)
  common <- names(table(c(repl1_cuts, repl2_cuts)))[table(c(repl1_cuts, repl2_cuts)) == 2]
  print(length(repl1_cuts))
  print(length(repl2_cuts))
  print(length(common))

  repl1 <- repl1 %>%
    filter(V1 %in% common)
  repl2 <- repl2 %>%
    filter(V1 %in% common)

  print("piRNAs making the cuts in replicate 1, replicate 2 and common for both replicates/samples (after filtering for common cuts):")
  repl1_pirnas <- unique(repl1$V4)
  repl2_pirnas <- unique(repl2$V4)
  common <- names(table(c(repl1_pirnas, repl2_pirnas)))[table(c(repl1_pirnas, repl2_pirnas)) == 2]
  print(length(repl1_pirnas))
  print(length(repl2_pirnas))
  print(length(common))

  tab <- rbind(repl1, repl2)
} else {
  tab <- repl1
}

if (!secondary) {
  tab <- tab %>%
    select(-V5) %>%
    distinct() # Remove expressions
} else {
  tab <- tab %>%
    distinct() # Remove expression
}

rio::export(x = tab, file = otab, format = "tsv", col.names = F)