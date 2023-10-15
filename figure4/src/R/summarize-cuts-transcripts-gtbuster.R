#!/usr/bin/env Rscript
#
# Visualize positions of predicted cuts in piRNA target cutting prediction in overall transcripts
#
################################################################################

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("plyranges"))
library("GenomicRanges")
suppressPackageStartupMessages(library("dplyr"))
library("ggplot2")

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

# Make binnning function
add_bin <- function(start, end, position) {
  bin <- cut(
    x = seq(start + 10, end + 10, by = 1),
    breaks = 20,
    labels = seq(1, 20, by = 1), include.lowest = T
  )[position - start + 1]
  return(bin)
}

annotate_cuts <- function(tab) {
  # tab <- ifile_predict

  tab_granges <- tab %>%
    bed_to_granges()

  tab_annot <- annot_granges(granges = bed_to_granges(tab), annot = bed_to_granges(bed), annot_name = "region", ignore.strand = T, annot_field = "name")
  tab_annot <- unique(tab_annot)
  tab_annot <- tab_annot[!is.na(mcols(tab_annot)$region)] # Remove unannotated (non-coding, no CDS, ...). CHECK BEFORE REMOVING!

  tab_annot$transcript_id <- unlist(lapply(tab_annot$region, function(x) strsplit(x, split = ":", fixed = T)[[1]][1]))
  tab_annot$type <- unlist(lapply(tab_annot$region, function(x) strsplit(x, split = ":", fixed = T)[[1]][2]))
  tab_annot$type <- factor(tab_annot$type, levels = c("utr5", "cds", "utr3"))
  # Get summary per type
  print("Total number of transcripts:")
  print(length(unique(chrom(tab_annot))))

  return(tab_annot)
}

summarize_cuts <- function(tab_annot, annot_bed, ylim = c(0, 150)) {
  summary_types <- tab_annot %>%
    as.data.frame() %>%
    group_by(sample, type) %>%
    tally()
  print("Summary by types:")
  print(summary_types)

  summary_transcripts <- tab_annot %>%
    as.data.frame() %>%
    group_by(sample, transcript_id, type) %>%
    tally()
  print("Summary by transcripts:")
  print(summary_transcripts)

  # Merge with annot to get region sizes
  summary_transcripts$name <- paste(summary_transcripts$transcript_id, summary_transcripts$type, sep = ":")

  summary_transcripts <- summary_transcripts %>%
    left_join(as.data.frame(annot_bed))

  summary_transcripts$size <- summary_transcripts$end - summary_transcripts$start

  # Make heatmaps https://themockup.blog/posts/2020-08-28-heatmaps-in-ggplot2/
  # It's not easy to make independent fill gradients in facets; here is a workaround https://stackoverflow.com/questions/33907998/assigning-individual-high-and-low-fill-values-using-geom-tile-facet-wrap?rq=1
  hex_plot_utr5 <- lapply(unique(summary_transcripts$sample), function(x) {
    hex_plot <- ggplot(data = subset(summary_transcripts, sample == x & type == "utr5"), aes(x = size, y = n)) +
      geom_hex(
        binwidth = c(6000 / (max(ylim) * 2), 1)
      ) +
      scale_fill_gradient(low = "red", high = "yellow") +
      coord_cartesian(xlim = c(0, 6000), ylim = ylim) +
      ggtitle(paste0(x, " - 5 UTR")) +
      xlab("Size of the feature. Limited to 6000 nt. Binned.") +
      ylab("Number of cut positions.") +
      theme_bw()
  })

  hex_plot_cds <- lapply(unique(summary_transcripts$sample), function(x) {
    hex_plot <- ggplot(data = subset(summary_transcripts, sample == x & type == "cds"), aes(x = size, y = n)) +
      geom_hex(
        # geom_bin2d(
        binwidth = c(15000 / (max(ylim) * 2), 1)
      ) +
      scale_fill_gradient(low = "red", high = "yellow") +
      coord_cartesian(xlim = c(0, 15000), ylim = ylim) +
      ggtitle(paste0(x, " - CDS")) +
      xlab("Size of the feature. Limited to 15000 nt. Binned.") +
      ylab("Number of cut positions.") +
      theme_bw()
  })

  hex_plot_utr3 <- lapply(unique(summary_transcripts$sample), function(x) {
    hex_plot <- ggplot(data = subset(summary_transcripts, sample == x & type == "utr3"), aes(x = size, y = n)) +
      geom_hex(
        # geom_bin2d(
        binwidth = c(15000 / (max(ylim) * 2), 1)
      ) +
      scale_fill_gradient(low = "red", high = "yellow") +
      coord_cartesian(xlim = c(0, 15000), ylim = ylim) +
      ggtitle(paste0(x, " - 3 UTR")) +
      xlab("Size of the feature. Limited to 15000 nt. Binned.") +
      ylab("Number of cut positions.") +
      theme_bw()
  })

  do.call(gridExtra::grid.arrange, c(hex_plot_utr5, ncol = length(unique(summary_transcripts$sample))))
  do.call(gridExtra::grid.arrange, c(hex_plot_cds, ncol = length(unique(summary_transcripts$sample))))
  do.call(gridExtra::grid.arrange, c(hex_plot_utr3, ncol = length(unique(summary_transcripts$sample))))

  return(summary_transcripts)
}

################################################################################

option_list <- list(
  make_option(
    c("-i", "--predict"),
    type = "character",
    help = "Predicted cuts only from piRNA sequence, not considering degradome-seq.", metavar = "File"
  ),
  make_option(
    c("-l", "--real"),
    type = "character",
    help = "Real cuts from piRNA sequence and degradome-seq data.", metavar = "File"
  ),
  make_option(
    c("-a", "--pirna"),
    type = "character",
    help = "File with piRNAs we want to subset; in our case, it's only piRNAs present in both replicates (either Het or RK).", metavar = "File"
  ),
  make_option(
    c("-b", "--bed"),
    type = "character",
    help = "Merged 5utr, cds, and 3utr annotation in BED format.", metavar = "File"
  ),
  make_option(
    c("-g", "--genetotrans"),
    type = "character",
    help = "Gene to transcript table for converting DE genes.", metavar = "File"
  ),
  make_option(
    c("-e", "--expr"),
    type = "character",
    help = "List of all expressed genes for cleaning unusued transcripts; we require at least 100 reads to map to the transcript to be considered expressed.", metavar = "File"
  ),
  make_option(
    c("-t", "--target_theor"),
    type = "character",
    help = "List of genes with predicted cuts but without degradome support.", metavar = "File"
  ),
  make_option(
    c("-n", "--trans_target"),
    type = "character",
    help = "List of transcripts with predicted cuts (with degradome/piRNA support) but not filtered at all", metavar = "File"
  ),
  make_option(
    c("-d", "--trans_de"),
    type = "character",
    help = "Table of differentially expressed genes.", metavar = "File"
  ),
  make_option(
    c("-r", "--de_target"),
    type = "character",
    help = "Selected differentially expressed genes which are targeted by piRNA/degradome combination but not filtered for RNA up/piRNA down yet and might not have strict filtering.", metavar = "File"
  ),
  make_option(
    c("-p", "--pval_mrna"),
    type = "double", default = 0.1,
    help = "Adjust pvalue for mRNA DE filtering; put 1.0 for no adj.p-value filtering. Default: 0.1"
  ),
  make_option(
    c("-q", "--pval_pirna"),
    type = "double", default = 0.25,
    help = "Adjust pvalue for piRNA DE filtering; put 1.0 for no adj.p-value filtering. Previously, we used 0.25 for the \"strict\" version. Default: 0.25"
  ),
  make_option(
    c("-o", "--odir"),
    type = "character",
    help = "Output dir for all the figures.", metavar = "Dir"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

ifile_predict <- opt$predict
ifile_real <- opt$real
subset_pirna <- opt$pirna
bed_merged <- opt$bed
gene_to_transcript <- opt$genetotrans
expressed_all <- opt$expr
trans_theo_target <- opt$target_theor
trans_target <- opt$trans_target
trans_de <- opt$trans_de
trans_de_target_regul <- opt$de_target
pval_mrna <- opt$pval_mrna
pval_pirna <- opt$pval_pirna
odir <- opt$odir

################################################################################
bed <- rio::import(bed_merged, format = "tsv")

bed <- bed %>%
  dplyr::rename("chr" = V1, "start" = V2, "end" = V3, "name" = V4, "score" = V5, "strand" = V6)

ifile_predict <- rio::import(ifile_predict, format = "tsv", header = T)
ifile_real <- rio::import(ifile_real, format = "tsv", header = T) # Should be bed already

subset_pirna <- rio::import(subset_pirna, format = "tsv", header = F)

dir.create(odir, recursive = T)

ifile_predict <- ifile_predict %>%
  rename("chr" = V1, "end" = V2) %>%
  mutate("start" = end - 1, "name" = ".", "score" = "0", "strand" = "+", "pirna" = V4)

ifile_real <- ifile_real %>%
  mutate(score = 0)

# Keep only piRNAs in our subset
ifile_predict <- ifile_predict %>%
  filter(V4 %in% !!subset_pirna$V1)
ifile_real <- ifile_real %>%
  filter(pirna %in% !!subset_pirna$V1)

ifile_predict_annot <- annotate_cuts(ifile_predict)
pdf(paste(odir, "size_vs_cuts-predict.pdf", sep = "/"), width = 8, height = 4)
summary_transcripts_predict <- summarize_cuts(tab_annot = ifile_predict_annot, annot_bed = bed, ylim = c(0, 150))
dev.off()

ifile_real_annot <- annotate_cuts(ifile_real) # missing transcripts could be transcripts without an annotated 5UTR/CDS/3UTR

pdf(paste(odir, "size_vs_cuts-real.pdf", sep = "/"), width = 8, height = 4)
summary_transcripts_real <- summarize_cuts(tab_annot = ifile_real_annot, annot_bed = bed, ylim = c(0, 25))
dev.off()

################################################################################

plot_cut_distro <- function(tab, subset = NULL, title = "", binary = FALSE, plotme = TRUE, plot_only = FALSE) {
  if (!is.null(subset)) { # without subsetting
    print("Subsetting")
    tab <- tab %>%
      filter(transcript_id %in% !!subset)
  }

  if (!plot_only) {
    if (binary) {
      print("Using binary \"counts\"")
      tab$total_trans_type_bin[tab$total_trans_type_bin > 0] <- 1
      tab$n_sizeNorm <- tab$total_trans_type_bin # If it's binary, don't normalize to bin length
    } else {
      # Normalize counts to bin size
      tab$n_sizeNorm <- tab$total_trans_type_bin / tab$bin_size
    }

    tab <- tab %>%
      group_by(sample, type, bin) %>%
      mutate(total_type_bin_sizeNorm = sum(n_sizeNorm)) # Get total sum per bin per type
  }

  # Plot raw, size normalized counts
  p1 <- ggplot(tab %>% select(sample, type, bin, total_type_bin_sizeNorm) %>% distinct(), aes(x = bin, y = total_type_bin_sizeNorm, color = sample, group = sample)) + # , shape=type
    geom_point() +
    scale_shape_manual(values = c(0, 1, 2)) +
    geom_line() +
    facet_grid(~type) +
    theme_classic() +
    theme(legend.position = "bottom") +
    ggtitle(paste0(title, "Absolute counts normalized to bin length."))

  ## Plot relative to number of cuts and to number of transcripts
  # Make total_type counts relative to number of transcripts - easier when comparing multiple sets
  if (!plot_only) {
    tab <- tab %>%
      group_by(sample) %>%
      mutate(total_sample = sum(total_trans_type_bin)) %>%
      group_by(sample, type, bin) %>%
      mutate(total_sizeNorm_rel = total_type_bin_sizeNorm / total_sample) %>% # sort of RPKM normalization; has to be to TOTAL number of reads per sample, otherwise if done per transcript it penalizes transcripts with a lot of cuts
      ungroup()
    tab <- tab %>%
      mutate(total_type_rel = total_sizeNorm_rel / length(unique(transcript_id))) # normalize to number of transcripts for all samples together
  }

  p2 <- ggplot(tab %>% select(sample, type, bin, total_type_rel) %>% distinct(), aes(x = bin, y = total_type_rel, color = sample, group = sample)) + # , shape=type
    geom_point() +
    scale_shape_manual(values = c(0, 1, 2)) +
    geom_line() +
    facet_grid(~type) +
    theme_classic() +
    theme(legend.position = "bottom") +
    ggtitle(paste0(title, "Relative counts (absolute counts normalized to bin length\n& total number of transcript cut positions & number of transcripts)."))

  if (plotme) {
    print(p1)
    print(p2)
  }

  return(tab)
}

# Plot distribution along the features
plot_distro <- function(ifile, ifile_summary, binary = FALSE, subset = NULL, plotme = TRUE, plot_only = FALSE) {
  if (!is.null(subset)) {
    print("Subsetting")
    ifile <- ifile %>%
      filter(transcript_id %in% !!subset)
    ifile_summary <- ifile_summary %>%
      filter(transcript_id %in% !!subset)
  }

  if (!plot_only) {
    ifile_annot <- ifile %>%
      as.data.frame() %>%
      rename("position" = end) %>%
      select(sample, transcript_id, position, type, pirna)
    ifile_annot <- ifile_annot %>%
      left_join(ifile_summary) %>%
      filter(size >= 60) # At least 60 nt in size (20 bins by 3 nt)

    # Add bins
    ifile_annot$bin <- apply(ifile_annot[, c("start", "end", "position")], 1, function(x) add_bin(start = x[1], end = x[2], position = x[3]))
    ifile_annot$bin_size <- ceiling(ifile_annot$size / 20) # Approx bin size for normalization

    # Make summary relative to total number of cuts per transcript
    coverage_features <- ifile_annot %>%
      group_by(sample, transcript_id) %>%
      mutate(total = n()) %>% # Get total number of cuts per transcript and normalize to norm_factor
      group_by(sample, transcript_id, type, bin) %>%
      mutate(total_trans_type_bin = n()) %>% # Get number of cuts per bin per transcript-type
      mutate(perc = total_trans_type_bin / total) %>% # Get percentage of cuts per bin
      select(sample, transcript_id, type, bin, bin_size, total, total_trans_type_bin, perc) %>%
      distinct() %>%
      ungroup() # Summarize

    coverage_features$type <- factor(coverage_features$type, levels = c("utr5", "cds", "utr3"))
  } else {
    coverage_features <- ifile
  }

  distro <- plot_cut_distro(coverage_features, subset = subset, binary = binary, plotme = plotme, plot_only = plot_only)
  return(distro)
}

#################################################################################

gene_to_transcript <- read.table(gene_to_transcript, header = T, stringsAsFactors = F, sep = "\t", quote = "\"") %>%
  select(gene_id, transcript_id)

expressed_all_list <- read.table(expressed_all, header = F, stringsAsFactors = F) %>%
  pull(V1) %>%
  unique() # All expressed transcripts
trans_theo_targeted_list <- read.table(trans_theo_target, header = F, stringsAsFactors = F) %>%
  pull(V1) %>%
  unique() # All theoretically targeted transcripts (no degradome information)
trans_target_list <- read.table(trans_target, header = T, stringsAsFactors = F, sep = "\t", quote = "\"") %>%
  left_join(gene_to_transcript) %>% # Convert to transcripts
  filter(transcript_id %in% expressed_all_list) %>% # Get only known expressed transcripts to avoid extras in the gene->transcript conversion
  pull(transcript_id) %>%
  unique() # All targeted but w/o information from DE (piRNA/RNA-Seq)
trans_de_list <- read.table(trans_de, header = T, stringsAsFactors = F, sep = "\t", quote = "\"") %>%
  filter(abs(logFC) >= log2(2) & PValueAdj < !!pval_mrna) %>% # Get transcripts
  left_join(gene_to_transcript) %>% # Convert to transcripts
  filter(transcript_id %in% expressed_all_list) %>% # Get only known expressed transcripts to avoid extras in the gene->transcript conversion
  pull(transcript_id) %>%
  unique() # All differentially expressed transcripts
trans_notde_target_list <- read.table(trans_target, header = T, stringsAsFactors = F, sep = "\t", quote = "\"") %>%
  filter(PValueAdj >= !!pval_mrna) %>%
  left_join(gene_to_transcript) %>% # Convert to transcripts
  filter(transcript_id %in% expressed_all_list) %>% # Get only known expressed transcripts to avoid extras in the gene->transcript conversion
  pull(transcript_id) %>%
  unique() # All NOT differentially expressed but targeted (but w/o information from piRNA)
trans_de_target_list <- read.table(trans_target, header = T, stringsAsFactors = F, sep = "\t", quote = "\"") %>%
  filter(PValueAdj < !!pval_mrna) %>%
  left_join(gene_to_transcript) %>% # Convert to transcripts
  filter(transcript_id %in% expressed_all_list) %>% # Get only known expressed transcripts to avoid extras in the gene->transcript conversion
  pull(transcript_id) %>%
  unique() # All differentially expressed but targeted but w/o information from piRNA
trans_de_target_regul_list <- read.table(trans_de_target_regul, header = T, stringsAsFactors = F, sep = "\t", quote = "\"") %>%
  filter(PValueAdj < !!pval_mrna & PValueAdj_pirna < !!pval_pirna) %>%
  filter(change_pirna == "Down" & change_rna == "Up") %>% # Get genes
  left_join(gene_to_transcript) %>% # Convert to transcripts
  filter(transcript_id %in% expressed_all_list) %>% # Get only known expressed transcripts to avoid extras in the gene->transcript conversion
  pull(transcript_id) %>%
  unique() # All differentially expressed and piRNA regulated and in correct orientation
trans_notde_target_regul_list <- read.table(trans_target, header = T, stringsAsFactors = F, sep = "\t", quote = "\"") %>%
  filter(PValueAdj >= !!pval_mrna & PValueAdj_pirna < !!pval_pirna) %>%
  filter(change_pirna == "Down") %>% # & change_rna == "Up") %>%  # Get genes
  left_join(gene_to_transcript) %>% # Convert to transcripts
  filter(transcript_id %in% expressed_all_list) %>% # Get only known expressed transcripts to avoid extras in the gene->transcript conversion
  pull(transcript_id) %>%
  unique() # All NOT differentially expressed but piRNA regulated and piRNA goes down

### Annotate transcripts
transcripts_desc <- as.data.frame(expressed_all_list)
transcripts_desc$target_theo <- transcripts_desc$expressed_all_list %in% trans_theo_targeted_list
transcripts_desc$de <- transcripts_desc$expressed_all_list %in% trans_de_list
transcripts_desc$target <- transcripts_desc$expressed_all_list %in% trans_target_list
transcripts_desc$target_regul <- transcripts_desc$expressed_all_list %in% trans_de_target_regul_list
transcripts_desc$target_false_regul <- transcripts_desc$expressed_all_list %in% trans_notde_target_regul_list
colSums(transcripts_desc[, 2:ncol(transcripts_desc)])

### Analysis separately for
# a) expressed theoretically targeted without de or degradome/pirna target
# b) de without degradome/pirna target
# c) degradome/pirna target without de
# d) degradome/pirna target de without piRNA regulated
# e) degradome/pirna target de piRNA regulated

###
# Theoretically targeted genes
length(unique(gene_to_transcript[gene_to_transcript$transcript_id %in% unique(summary_transcripts_predict$transcript_id), "gene_id"]))
# Theoretically targeting piRNAs
length(unique(ifile_predict_annot$pirna))
# Degradome-seq targeted genes
length(unique(gene_to_transcript[gene_to_transcript$transcript_id %in% unique(summary_transcripts_real$transcript_id), "gene_id"]))
# Degradome-seq targeting piRNAs
length(unique(ifile_real_annot$pirna))
# Actually targeted genes
read.table(trans_de_target_regul, header = T, stringsAsFactors = F, sep = "\t", quote = "\"") %>%
  filter(PValueAdj < !!pval_mrna & PValueAdj_pirna < !!pval_pirna) %>%
  filter(change_pirna == "Down" & change_rna == "Up") %>% # Get genes
  left_join(gene_to_transcript) %>% # Convert to transcripts
  pull(gene_id) %>%
  unique() %>%
  length()
# Actually-seq targeting piRNAs
read.table(trans_de_target_regul, header = T, stringsAsFactors = F, sep = "\t", quote = "\"") %>%
  filter(PValueAdj < !!pval_mrna & PValueAdj_pirna < !!pval_pirna) %>%
  filter(change_pirna == "Down" & change_rna == "Up") %>% # Get genes
  left_join(gene_to_transcript) %>% # Convert to transcripts
  filter(transcript_id %in% expressed_all_list) %>%
  pull(pirna_merged) %>%
  strsplit(x = ., split = "|", fixed = T) %>%
  unlist() %>%
  unique() %>%
  length()

### Basic summary of how many transcripts are theoretically targeted from all the expressed
print("Total number of expressed transcripts")
nrow(transcripts_desc)
print("Total number of theoretically targeted transcripts and the percentage of total")
sum(transcripts_desc$target_theo)
round(sum(transcripts_desc$target_theo) / (nrow(transcripts_desc) / 100), 2)
print("Number of targeted transcripts (degradome+piRNA but unknown mRNA DE regulation and the percentage of total")
sum(transcripts_desc$target)
round(sum(transcripts_desc$target) / (nrow(transcripts_desc) / 100), 2)
print("Number of targeted transcripts (degradome+piRNA+mRNA DE regulation and the percentage of total")
sum(transcripts_desc$target_regul)
round(sum(transcripts_desc$target_regul) / (nrow(transcripts_desc) / 100), 2)

ifile_predict_annot_plot <- ifile_predict_annot %>%
  as.data.frame() %>%
  select(sample, transcript_id, start, end, type, pirna)
ifile_real_annot_plot <- ifile_real_annot %>%
  as.data.frame() %>%
  select(sample, transcript_id, start, end, type, pirna)

## Plot all transcripts together -  all theoretical and all predicted
# "Absolute"
plot_distro_predict <- plot_distro(ifile_predict_annot_plot, summary_transcripts_predict, plotme = FALSE) # predicted theoretical cuts

pdf(paste(odir, "predicted_cuts.pdf", sep = "/"), width = 8)
  plot_distro(plot_distro_predict, plot_only = T) # predicted theoretical cuts
dev.off()

plot_distro_real <- plot_distro(ifile_real_annot_plot, summary_transcripts_real, plotme = FALSE) # real cuts

pdf(paste(odir, "real_cuts.pdf", sep = "/"), width = 8)
  plot_distro(plot_distro_real, plot_only = T) # predicted theoretical cuts
dev.off()

pdf(paste(odir, "predicted_and_real_cuts.pdf", sep = "/"), width = 8)
  plot_distro_both <- plot_distro(ifile = rbind(plot_distro_predict, plot_distro_real), plot_only = T) # both real and predicted cuts
dev.off()

# Binary only - cut yes/no, no numbers
plot_distro_predict_bin <- plot_distro(ifile_predict_annot_plot, summary_transcripts_predict, binary = T, plotme = FALSE) # predicted theoretical cuts

pdf(paste(odir, "predicted_cuts-binary.pdf", sep = "/"), width = 8)
  plot_distro(plot_distro_predict_bin, plot_only = T) # predicted theoretical cuts
dev.off()

plot_distro_real_bin <- plot_distro(ifile_real_annot_plot, summary_transcripts_real, binary = T, plotme = FALSE) # real cuts

pdf(paste(odir, "real_cuts-binary.pdf", sep = "/"), width = 8)
  plot_distro(plot_distro_real_bin, plot_only = T) # predicted theoretical cuts
dev.off()

pdf(paste(odir, "predicted_and_real_cuts-binary.pdf", sep = "/"), width = 8)
  plot_distro_both_bin <- plot_distro(ifile = rbind(plot_distro_predict_bin, plot_distro_real_bin), plot_only = T) # both real and predicted cuts
dev.off()