#!/usr/bin/env Rscript
#
# Merge annotated GTBuster with results from DE expression - genes and piRNAs
#
# selected.tsv output - We filter by RNA diff. expression (adj.pval) but not piRNA diff. expression! Unless you change pval_pirna from the default 1.0
#
################################################################################

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ggrepel"))
library("ggplot2")
suppressPackageStartupMessages(library("dplyr"))

################################################################################

strip_gtbuster <- function(tab) {
  tab <- tab %>%
    select(-seqnames, -start, -end, -width, -strand, -pirna_pos)
  return(tab)
}

###############################################################################
option_list <- list(
  make_option(
    c("-g", "--gene_de"),
    type = "character",
    help = "RNA-Seq DE table for genes.", metavar = "File"
  ),
  make_option(
    c("-a", "--adjpval_rna"),
    type = "double", default = 0.1,
    help = "Adjusted p-value for mRNA DE filtering. Default: 0.1", metavar = "float",
  ),
  make_option(
    c("-e", "--gtbuster_het"),
    type = "character",
    help = "Filtered and annotated GTBuster cuts from the Het condition.", metavar = "File"
  ),
  make_option(
    c("-r", "--gtbuster_rk"),
    type = "character",
    help = "Filtered and annotated GTBuster cuts from the RK condition.", metavar = "File"
  ),
  make_option(
    c("-n", "--pirna_de"),
    type = "character",
    help = "Predicted targeting piRNA DE table.", metavar = "File"
  ),
  make_option(
    c("-p", "--adjpval_pirna"),
    type = "double", default = 1.0,
    help = "Adjusted p-value for piRNA DE filtering. Default: 1.0", metavar = "float",
  ),
  make_option(
    c("-d", "--degra_cuts"),
    type = "character",
    help = "Degradome cuts for all the samples. Use ONLY those containing piRNA targeted cuts, not overall degradome cuts (they'll have mostly Xrn1 and other stuff)", metavar = "File"
  ),
  make_option(
    c("-o", "--ofile"),
    type = "character",
    help = "Output file", metavar = "File"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

gene_de <- opt$gene_de
gtbuster_het <- opt$gtbuster_het
gtbuster_rk <- opt$gtbuster_rk
pirna_de <- opt$pirna_de
degra_cuts <- opt$degra_cuts
ofile <- opt$ofile
pval <- opt$adjpval_rna # 0.1 # adj.p-value for mRNA filtering
pval_pirna <- opt$adjpval_pirna # 1.0 # adj.p-value for piRNA filtering; set to 1.0 if you DON'T want to filter by pirna adj.p-value; should we use 0.25??
min_abundancebin <- 1 # min abundance bin for mRNA

################################################################################

gene_de <- rio::import(gene_de, format = "tsv")
gtbuster_het <- rio::import(gtbuster_het, format = "tsv") %>%
  strip_gtbuster()
gtbuster_rk <- rio::import(gtbuster_rk, format = "tsv") %>%
  strip_gtbuster()
pirna_de <- rio::import(pirna_de, format = "tsv")
degra_cuts <- rio::import(degra_cuts, format = "tsv")

gtbuster_het <- gtbuster_het %>%
  rename(pirna_mols_Het = pirna_mols, pirna_seeds_Het = pirna_seeds)
gtbuster_rk <- gtbuster_rk %>%
  rename(pirna_mols_RK = pirna_mols, pirna_seeds_RK = pirna_seeds)

gtbuster <- full_join(gtbuster_het, gtbuster_rk)

tab <- gene_de %>%
  left_join(gtbuster) %>%
  distinct()

pirna_de <- pirna_de %>%
  select(-change) %>%
  rename("PValueAdj" = FDR) %>%
  rename_all(., list(~ paste(., "pirna", sep = "_"))) %>%
  rename(pirna_merged = name_pirna, pirna_seed = pirna_seed_pirna, gene_id = gene_id_pirna)

print("foo")
tab <- tab %>%
  left_join(pirna_de) %>%
  select(-pirna, -LR, -PValue, -F_pirna, -PValue_pirna) %>%
  distinct()

# We might have some duplicated rows. This will happen if the original GTBuster file contained unique piRNAs for individual samples but overall the counts are the same.
# Check which piRNAs are duplicated and select rows with no NA values
# The duplication must be across degradation cut + piRNA merged name, not just name (we might have, in rare occasions, same combination of piRNAs at the same gene at different positions)
tab$helper <- paste(tab$degra_coord, tab$pirna_merged, sep = ",") # Make helper column for duplication

pirna_names <- tab$helper[duplicated(tab$helper)]
pirna_names <- pirna_names[!(pirna_names %in% "NA,NA")] # Remove NA,NA = no matching degra+piRNA

# Note: These are piRNAs common for all samples
tab_uniq <- tab %>%
  filter(!(helper %in% pirna_names))

# Get duplicated rows
tab_dupl_main <- tab %>%
  filter(helper %in% pirna_names)
# Get rows which don't have any NA values from the duplicated rows
# IMPORTANT: Hardcoded part of the following  - assuming two conditions in the input
tab_dupl_main$nas <- rowSums(is.na(tab_dupl_main[c("pirna_mols_Het", "pirna_seeds_Het", "pirna_mols_RK", "pirna_seeds_RK")]))

tab_dupl <- tab_dupl_main %>%
  filter(nas == 0) %>%
  select(-nas)

# Get the rest of duplicated rows
tab_dupl2 <- tab_dupl_main %>%
  filter(!(helper %in% tab_dupl$helper))

# If there was at least one common piRNA the rows should be already included in one of the previous outputs
# If there are not, we can sum the rest per degra cut and piRNA name
tab_dupl2 <- tab_dupl2 %>%
  group_by(degra_coord, pirna_merged) %>%
  mutate(
    pirna_mols_Het = sum(pirna_mols_Het, na.rm = T),
    pirna_seeds_Het = sum(pirna_seeds_Het, na.rm = T),
    pirna_mols_RK = sum(pirna_mols_RK, na.rm = T),
    pirna_seeds_RK = sum(pirna_seeds_RK, na.rm = T),
  ) %>%
  select(-nas) %>%
  distinct() %>%
  as.data.frame()

# Merge them back to one table
tab_out <- rbind(tab_uniq, tab_dupl)
tab_out <- rbind(tab_out, tab_dupl2)
tab_out <- tab_out %>%
  select(-helper)

# Get only targeted genes
tab_target <- tab_out %>%
  filter(!(is.na(pirna_merged)))

# Note: I checked, and duplicated genes/piRNA names are only in case we have multiple degradome target site in the same gene

# Potential target genes

### Visualization
# Add change based on logFC and very mild adj-pvalue for piRNA (because our piRNA replicates suck) for visualization
tab_target$change_pirna <- ""
tab_target$change_pirna[tab_target$logFC_pirna > 0 & tab_target$PValueAdj_pirna < pval_pirna] <- "Up"
tab_target$change_pirna[tab_target$logFC_pirna < 0 & tab_target$PValueAdj_pirna < pval_pirna] <- "Down"
tab_target$change_pirna[is.na(tab_target$change_pirna)] <- "Not DE"
tab_target$change_pirna <- factor(tab_target$change_pirna, levels = c("Down", "Up", "Not DE"))

# All results

# Only DE based on p-value
# Note: If there are genes more than once it means they have more target sites
p1 <- ggplot(
  subset(tab_target, PValueAdj < pval & AbundanceBin >= min_abundancebin & AbundanceBin_pirna >= 1),
  aes(x = logFC, y = logFC_pirna, color = change_pirna, shape = type)) +
  geom_point(alpha = 0.7) + # shape=16,
  geom_hline(
    yintercept = 0, linetype = "dashed",
    color = "black", size = 0.5
  ) +
  geom_vline(
    xintercept = 0, linetype = "dashed",
    color = "black", size = 0.5
  ) +
  theme_bw() +
  geom_text_repel(
    data = subset(tab_target, PValueAdj < pval & AbundanceBin >= min_abundancebin & AbundanceBin_pirna >= 1),
    aes(label = gene_name),
    size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )

dir.create(dirname(ofile), recursive = T)

pdf(file = gsub(".tsv", ".pdf", ofile), width = 8)
  print(p1)
dev.off()

# We filter by RNA diff. expression (adj.pval) but not piRNA diff. expression (except very unlikely piRNAs to be DE)!
selected <- subset(tab_target, PValueAdj < pval & AbundanceBin >= min_abundancebin & PValueAdj_pirna < pval_pirna & AbundanceBin_pirna >= 1)
selected$change_rna <- NA
selected$change_rna[selected$logFC > 0] <- "Up"
selected$change_rna[selected$logFC < 0] <- "Down"

# Add degradome cut support reads
degra_cuts <- degra_cuts %>% # Get degradome cuts RPM, not normalized to mRNA abundance!!!
  tidyr::pivot_longer(cols = starts_with("mmu."), names_to = "sample", values_to = "count") %>%
  group_by(sample) %>%
  mutate(rpm = round(count / (sum(count, na.rm = T) / 1000000), 3)) %>%
  tidyr::pivot_wider(names_from = sample, values_from = c(count, rpm)) %>%
  select(name, starts_with("rpm_"))

selected <- selected %>%
  left_join(
    degra_cuts %>%
      rename("degra_coord" = name) %>%
      rename_all( # Rename Degradome-Seq
        list(~
          gsub("mmu.PARESeq.polya.", "degra_", .))
      ),
    by = "degra_coord"
  ) %>%
  rename_all( # rename RNA-Seq
    list(~
      gsub("mmu.RNASeq.total.", "", .))
  )


selected <- selected %>%
  arrange(change_pirna, change_rna)

# Get degraome occupancy normalized by mRNA levels
selected$degra_abund_normByRNA_MiwiHet.P24 <-
  rowMeans(selected[, c(grep("^rpm_degra_MiwiHet.P24.", colnames(selected)))]) /
    rowMeans(selected[, c(grep("^MiwiHet.P24.*-abundance$", colnames(selected)))])
selected$degra_abund_normByRNA_MiwiRK.P24 <-
  rowMeans(selected[, c(grep("^rpm_degra_MiwiRK.P24.", colnames(selected)))]) /
    rowMeans(selected[, c(grep("^MiwiRK.P24.*-abundance$", colnames(selected)))])

selected$logFC_degra_normByRNA <- log2(selected$degra_abund_normByRNA_MiwiRK.P24 / selected$degra_abund_normByRNA_MiwiHet.P24)

write.table(x = selected, file = gsub(".tsv", "-selected.tsv", ofile), sep = "\t", row.names = F, quote = F)
write.table(x = tab_target, file = ofile, sep = "\t", row.names = F, quote = F)