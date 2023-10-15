#!/usr/bin/env Rscript
#
# Conversion of genomic to transcriptomic coordinates
#
# main source: https://www.biostars.org/p/251076/
# for data.frame: https://www.biostars.org/p/230758/
#
################################################################################

library("optparse")
suppressPackageStartupMessages(library("GenomicFeatures")) # For general gtf, not ensembl

# Read command line options and arguments
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "Gene annotation in GTF to use for the conversion.", metavar = "File"
  ),
  make_option(
    c("-o", "--output"),
    type = "character", default = "stdout",
    help = "Genomic elements relative to the transcript. Will add :mrna to transcripts with cds and :ncrna without.
      Will annotate :cds, :utr3, :utr5 where available.Can be stdout (use 'stdout').", metavar = "File"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

##################################################################################
txdb <- makeTxDbFromGFF(file = opt$input, format = "gtf") # load db
tx_lens <- transcriptLengths(txdb, with.cds_len = TRUE, with.utr5_len = TRUE, with.utr3_len = TRUE) # get length of transcript, cds, 3utr, 5utr

# Remove tx where tx length is not equal to cds+utr5+utr3 length where all of them are present so we don't get dual annotation
index <- tx_lens$cds_len > 0 & tx_lens$utr5_len > 0 & tx_lens$utr3_len > 0
tx_lens_com <- tx_lens[index, ]
tx_lens_incom <- tx_lens[!index, ]
index2 <- tx_lens_com$tx_len == tx_lens_com$cds_len + tx_lens_com$utr5_len + tx_lens_com$utr3_len
tx_lens_com <- tx_lens_com[index2, ]

tx_lens <- rbind(tx_lens_com, tx_lens_incom)

# UTRs and CDS
tx_lens$annot <- NULL
tx_lens$start <- 0
tx_lens$stop <- 0

tx_lens_cds <- tx_lens[tx_lens$cds_len > 0, ]
tx_lens_utr5 <- tx_lens[tx_lens$utr5_len > 0, ]
tx_lens_utr3 <- tx_lens[tx_lens$utr3_len > 0, ]

tx_lens_cds$start <- tx_lens_cds$utr5_len
tx_lens_cds$stop <- tx_lens_cds$utr5_len + tx_lens_cds$cds_len

tx_lens_utr5$stop <- tx_lens_utr5$utr5_len

tx_lens_utr3$start <- tx_lens_utr3$utr5_len + tx_lens_utr3$cds_len
tx_lens_utr3$stop <- tx_lens_utr3$start + tx_lens_utr3$utr3_len

tx_lens_cds$annot <- paste(tx_lens_cds$tx_name, "cds", sep = ":")
tx_lens_utr5$annot <- paste(tx_lens_utr5$tx_name, "utr5", sep = ":")
tx_lens_utr3$annot <- paste(tx_lens_utr3$tx_name, "utr3", sep = ":")

# mRNA and ncRNA
tx_lens$stop <- tx_lens$tx_len
tx_lens$annot[tx_lens$cds_len == 0] <- paste(tx_lens$tx_name[tx_lens$cds_len == 0], "ncrna", sep = ":")
tx_lens$annot[tx_lens$cds_len > 0] <- paste(tx_lens$tx_name[tx_lens$cds_len > 0], "mrna", sep = ":")

# Merge
tx_lens_out <- rbind(tx_lens, tx_lens_cds, tx_lens_utr5, tx_lens_utr3)

# Add strand
# tx_lens_out<-merge(tx_lens_out, tx_full[c("tx_name", "strand")])

# Organize
# ENST00000335137	0	1054	ENST00000335137:mrna	1054	+ # original manolis structure
# tx_lens_out <- tx_lens_out[, c("tx_name", "start", "stop", "annot", "tx_len", "strand")]
tx_lens_out <- tx_lens_out[, c("tx_name", "start", "stop", "annot", "tx_len")]
tx_lens_out$strand <- "+"

# Write
if (opt$output == "stdout") {
  write.table(x = tx_lens_out, file = "", quote = F, sep = "\t", col.names = F, row.names = F)
} else {
  write.table(x = tx_lens_out, file = opt$output, quote = F, sep = "\t", col.names = F, row.names = F)
}
