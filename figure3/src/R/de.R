#!/usr/bin/env Rscript
#
# Calculate differential gene expression and some basic stats
#

suppressPackageStartupMessages(library(dplyr)) # Loaded later on because of conflicts with biomaRt
library(ggplot2) # Loaded later on because of conflicts with biomaRt
suppressPackageStartupMessages(library(gplots))
library(RColorBrewer)
library(grid)
# library(biomaRt) # Loaded only if the input GTF doesn't exist and we need to make fresh gene description file
suppressPackageStartupMessages(library(rtracklayer))
library(tximport)
suppressPackageStartupMessages(library(edgeR))
# suppressMessages(library(tibble)) # Only used once to rename
# suppressMessages(library(reshape2)) # Only used once to melt df

################################################################################

args <- commandArgs(trailingOnly = TRUE)

# Allow multiple files at the same time
hh <- paste(unlist(args), collapse = " ")
listoptions <- unlist(strsplit(hh, "--"))[-1]
options.args <- sapply(listoptions, function(x) {
  unlist(strsplit(x, " "))[-1]
})
options.names <- sapply(listoptions, function(x) {
  option <- unlist(strsplit(x, " "))[1]
})

names(options.args) <- unlist(options.names)

print(options.args)

# Check number of arguments
if (length(options.args) < 1) {
  options.args <- character(1)
  names(options.args) <- "help"
} else if (length(options.args) < 11) {
  stop("Please specify all arguments.", call. = FALSE)
}

# Help section
if ("help" %in% names(options.args)) {
  cat("Usage:
    ./de.R --samples input file(s) --odir output plot in pdf ...

    For example: ./de.R \
    --samples mmu.RNASeq.Het.P24.1 mmu.RNASeq.Het.P24.2 mmu.RNASeq.Het.P24.3 mmu.RNASeq.RK.P24.1 mmu.RNASeq.RK.P24.2 mmu.RNASeq.RK.P24.3 \
    --idir /home/joppelt/projects/pirna_mouse/samples \
    --gtf /home/joppelt/projects/pirna_mouse/data/mm10/Mus_musculus.GRCm38.99.gtf \
    --odir /home/joppelt/projects/pirna_mouse/analysis/expression/results \
    --design /home/joppelt/projects/pirna_mouse/analysis/expression/data/design-rna.txt \
    --pval 0.1 \
    --fc 2 \
    --compcond1 Het \
    --compcond2 RK \
    --refcond Het \
    --counts salmon \
    --replicates TRUE \
    --batch TRUE
    ")
  q(save = "no")
}

################################################################################

bamlist <- options.args$samples
dir <- options.args$idir
ingtf <- options.args$gtf
resdir <- options.args$odir
design_tab <- options.args$design
pvalset <- as.numeric(options.args$pval)
fcset <- as.numeric(options.args$fc)
compcond1 <- options.args$compcond1 # compare condition 1
compcond2 <- options.args$compcond2 # compare condition 2
refcond <- options.args$refcond # reference condition
counts <- options.args$counts # [salmon|classic] using Salmon or "classic" counts such as STAR or featureCounts
ifile <- options.args$ifile # if we have 'classic' counts, we have should specific input file with merged counts; if NULL then we scan for *counts.merged.tab and hope it's there only once
batchCor <- options.args$batch # attempt batch correction if there is enough conditions/sample; main purpose is to disable batch correction if we don't want it

if (is.null(options.args$replicates)) {
  repl <- TRUE # Default if no value was given
} else {
  repl <- as.logical(options.args$replicates)
}
min_count <- 10 # so far used only for stats of total exp. genes

################################################################################
fcset <- log(fcset, 2)

# Default for counts
if (is.null(counts)) {
  counts <- "salmon"
}

resdir <- paste(resdir, "de", paste0(compcond2, "vs", compcond1), sep = "/")
dir.create(resdir, recursive = T)

coldata <- read.table(design_tab, sep = "\t", header = T, stringsAsFactors = F)

if (!file.exists(paste(resdir, "data", "mart-desc.tsv", sep = "/"))) {
  print("Existing gene description file was not found, trying to parse input GTF name and fetch description from biomaRt.")

  library(biomaRt)
  # This will work ok for human and mouse. For other organisms please change gene_name variable specific for that organism
  # Get Ensembl host
  ens <- biomaRt::listEnsemblArchives()
  ens_version <- strsplit(ingtf, ".", fixed = T)[[1]][3]
  ens_version <- gsub("-w_.*", "", ens_version) # Strip spikein from the version if it exists (we name ref. gtf with spikein as Mus_musculus.GRCm38.99-w_ercc.gtf, for example
  ens_host <- ens[ens$version == ens_version, "url"]
  ens_host <- gsub("http://", "", ens_host)

  # Get organism
  if (length(grep("Homo_sapiens", ingtf)) > 0) {
    organism <- "hsapiens_gene_ensembl"
    gene_name <- "hgnc_symbol"
  } else if (length(grep("Mus_musculus", ingtf)) > 0) {
    organism <- "mmusculus_gene_ensembl"
    gene_name <- "mgi_symbol"
  }

  m <- biomaRt::useMart(
    host = ens_host,
    biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = organism
  )

  # useCache = FALSE doesn't exist in some earlier versions of biomaRt package; I checked only 2.40.1 (doesn't have it) and 2.42.2 (has it)
  if (gsub(pattern = ".", replacement = "", packageVersion("biomaRt"), fixed = T) >= 2420) {
    mart_desc <- biomaRt::getBM(
      attributes = c("ensembl_gene_id", "ensembl_transcript_id", gene_name, "gene_biotype", "transcript_biotype", "description"), mart = m,
      useCache = FALSE
    ) # useCache = FALSE helps with filter_ error https://support.bioconductor.org/p/p132704/; some versions might issue an error with unused argument
  } else {
    mart_desc <- biomaRt::getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", gene_name, "gene_biotype", "transcript_biotype", "description"), mart = m) # useCache = FALSE helps with filter_ error https://support.bioconductor.org/p/p132704/; some versions might issue an error with unused argument
  }

  mart_desc <- mart_desc %>%
    dplyr::rename(gene_id = ensembl_gene_id, transcript_id = ensembl_transcript_id, gene_name = !!gene_name) # "entrezgene_id"; entrez_id = "entrezgene_id"
  mart_desc$gene_name[mart_desc$gene_name == ""] <- mart_desc$gene_id[mart_desc$gene_name == ""]
  mart_desc[, gene_name] <- NULL # Remove old column, if still exists
  mart_desc <- unique(mart_desc) # I think this used to be [, 1] pointing to some Ensembl numeric ID but it seems it's not there anymore

  if (!(file.exists(ingtf))) {
    print("Warning: Input GTF doesn't exist, we'll use biomaRt as gene description which might cause problems with Salmon tx2gene import as some transcripts might be missing.")
    dir.create(paste(resdir, "data", sep = "/"))
    write.table(mart_desc, paste(resdir, "data", "mart-desc.tsv", sep = "/"), sep = "\t", row.names = F)
  }
} else {
  print("Existing gene description file was found, will use that.")
}

# Remove unused samples
if (nrow(coldata) != length(bamlist)) {
  print("There are some unused samples (not present in the input list), will remove them from the design: ")
  coldata %>%
    filter(!sample %in% !!bamlist) %>%
    dplyr::select(sample) %>%
    unlist(use.names = F) %>%
    print()

  coldata <- coldata %>%
    filter(sample %in% !!bamlist)
}

if (!("condition" %in% colnames(coldata))) {
  print("No \"condition\" column in the input design table, guessing from the bamlist")
  coldata$condition <- bamlist %>%
    #    gsub("mmu.", "", .) %>%
    #    gsub("Seq.", "", .) %>%
    #    gsub("\\.([0-9])", "", ., fixed = F) %>%
    #    gsub("RNA|RIBO", "", .) %>%
    stringr::str_split(stringr::fixed("."), simplify = T) %>%
    as.data.frame() %>%
    dplyr::select(V3) %>%
    unlist()
} else {
  "Using \"condition\" column from input design table"
}

if (repl == FALSE) { # Set batches/patients
  print("\"repl\" set to false, filling up \"batch\" column randomly")
  coldata$batch <- c(1:nrow(coldata))
} else if (repl == TRUE) {
  if (!("batch" %in% colnames(coldata))) {
    print("No \"batch\" column in input design table, guessing from the bamlist")
    coldata$batch <- sub(".*\\.", "", bamlist)
  } else {
    "Using \"batch\" column from input design table"
  }
} else {
  stop("Please set --replicates to either TRUE or FALSE")
}

# make sure all are factors
coldata$condition <- factor(coldata$condition, levels = unique(coldata$condition))
coldata$batch <- factor(coldata$batch, levels = unique(coldata$batch))
coldata$condition <- relevel(coldata$condition, ref = refcond)
coldata$batch <- relevel(coldata$batch, ref = levels(coldata$batch)[1])

# Order bamlist to match coldata
bamlist <- bamlist[match(coldata$sample, bamlist)]
names(bamlist) <- coldata$name

print("Final design table check")
print(coldata)
write.table(x = coldata, file = paste(resdir, "design.txt", sep = "/"), sep = "\t", col.names = NA, quote = F)

if (counts == "salmon") {
  files <- file.path(dir, coldata$sample, "counts", "quant.sf")
  names(files) <- coldata$name
} else {
  if (is.null(ifile)) {
    files <- list.files(path = dir, pattern = "*counts.merged.tab", full.names = T)
    if (length(files) != 1) {
      stop("Cannot find input file or finding too many. Please use --ifile to specify the input file.")
    }
  } else {
    files <- ifile
  }
}

print("Input files check")
print(files)

if (!all(file.exists(files))) {
  stop("At least one of the count files doesn't exist, please check.")
}

###
## Analysis of differentially expressed genes with edgeR

# Get gene annotation (mainly to remove rRNA later on)
if (!file.exists(paste(resdir, "data", "mart-desc.tsv", sep = "/"))) {
  if (file.exists(ingtf)) {
    print("Reading and parsing input gtf for description file.")

    desc <- as.data.frame(rtracklayer::import(ingtf, format = "gtf")) # Load GTF

    desc$gene_name[is.na(desc$gene_name)] <- desc$gene_id[is.na(desc$gene_name)]

    desc <- desc %>%
      dplyr::select(gene_id, transcript_id, gene_name, gene_biotype, transcript_name, transcript_biotype) %>%
      filter(!is.na(transcript_id)) %>%
      distinct()

    mart_desc <- mart_desc %>%
      dplyr::select(-starts_with("gene_"), -transcript_biotype)
    desc <- desc %>%
      left_join(mart_desc)

    dir.create(paste(resdir, "data", sep = "/"))
    write.table(desc, paste(resdir, "data", "mart-desc.tsv", sep = "/"), sep = "\t", row.names = F)
  } else {
    stop("Cannot find input gtf and pre-generated gene annotation doesn't exist.")
  }
} else {
  print("Reading annotation file.")
  desc <- read.table(paste(resdir, "data", "mart-desc.tsv", sep = "/"), sep = "\t", stringsAsFactors = F, header = T)
}

tx2gene <- desc %>%
  dplyr::select(transcript_id, gene_id) %>%
  dplyr::rename(TXNAME = "transcript_id", GENEID = "gene_id") %>%
  distinct()
if (!file.exists(paste(resdir, "data", "tx2gene.tsv", sep = "/"))) {
  write.table(tx2gene, paste(resdir, "data", "tx2gene.tsv", sep = "/"), sep = "\t", row.names = F)
} else {
  tx2gene <- read.table(paste(resdir, "data", "tx2gene.tsv", sep = "/"), sep = "\t", stringsAsFactors = F, header = T)
}

# Get only gene annotation
desc <- desc %>%
  dplyr::select(-starts_with("transcript")) %>%
  distinct()

## Import counts
# Get rRNA genes for cleaning later
rrna <- desc %>%
  filter(grepl("*rRNA", gene_biotype)) %>%
  pull(gene_id)

if (counts == "salmon" | is.null(counts)) {
  # Import salmon quantification
  # QuantSeq recommendations are slightly different https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#3%E2%80%99_tagged_RNA-seq

  txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

  txi$abundance <- subset(txi$abundance, !(rownames(txi$abundance) %in% rrna))
  txi$counts <- subset(txi$counts, !(rownames(txi$counts) %in% rrna))
  txi$length <- subset(txi$length, !(rownames(txi$length) %in% rrna))

  # Remove 0 expressed genes
  keep <- rowSums(txi$counts) > 0
  txi$abundance <- txi$abundance[keep, ]
  txi$counts <- txi$counts[keep, ]
  txi$length <- txi$length[keep, ]
} else {
  # Import 'classic' counts (STAR/featureCounts/HTSeq)

  txi <- NULL
  txi$counts <- read.table(ifile, header = TRUE, row.names = 1)
  txi$counts <- txi$counts[!rownames(txi$counts) %in% c("N_ambiguous", "N_multimapping", "N_noFeature", "N_unmapped"), ] # Remove header - STAR counts
  txi$counts <- txi$counts[, !colnames(txi$counts) %in% c("Chr", "Start", "End", "Strand", "Length")] # # Remove header - featureCounts counts

  # Rename input files/columns
  txi$counts <- txi$counts[, coldata$sample] # reorder
  colnames(txi$counts) <- coldata$name # rename

  txi$counts <- subset(txi$counts, !(rownames(txi$counts) %in% rrna))

  # Remove 0 expressed genes
  keep <- rowSums(txi$counts) > 0
  txi$counts <- txi$counts[keep, ]
}

# Prepare counts for edgeR https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#edger
cts <- txi$counts
cts <- cts[, coldata$name] # Make sure cts and design are ordered the same

if (counts == "salmon" | is.null(counts)) {
  txi$abundance <- txi$abundance[, coldata$name] # Make sure abundance and design are ordered the same
  txi$length <- txi$length[, coldata$name] # Make sure length and design are ordered the same

  normMat <- txi$length

  # Obtaining per-observation scaling factors for length, adjusted to avoid
  # changing the magnitude of the counts.
  normMat <- normMat / exp(rowMeans(log(normMat)))
  normCts <- cts / normMat
  # Computing effective library sizes from scaled counts, to account for
  # composition biases between samples.
  eff.lib <- calcNormFactors(normCts) * colSums(normCts)
  # Combining effective library sizes with the length factors, and calculating
  # offsets for a log-link GLM.
  normMat <- sweep(normMat, 2, eff.lib, "*")
  normMat <- log(normMat)

  # Creating a DGEList object for use in edgeR.
  y <- DGEList(cts, group = coldata$condition) # genes=rownames(cts) would be more usefull if we would put here some external annotation and then filter according to that
  y <- scaleOffset(y, normMat)
} else {
  y <- DGEList(cts, group = coldata$condition) # genes=rownames(cts) would be more usefull if we would put here some
  txi$abundance <- cpm(y)[, colnames(cts)]
}

# filtering by "existence", annotation, and min. expression
isexpr <- rowSums(y$counts) > 0

keep <- filterByExpr(y, group = coldata$condition) # edgeR manual specifies we should not use "baseline" coldata for this filtering but only the "difference" condition
print("Number of genes before expression filtering:")
nrow(y)
y <- y[isexpr & keep, ]
print("Number of genes after expression filtering:")
nrow(y)
# y is now ready for estimate dispersion functions see edgeR User's Guide

# Get stats of # of expressed genes and exp. distribution
exp_genes <- reshape2::melt(as.matrix(cts)) %>%
  group_by(Var2) %>%
  filter(value > !!min_count) %>%
  summarise(exp_genes = n(), raw_expr = sum(value))

p <- ggplot(exp_genes, aes(x = Var2, y = exp_genes, fill = Var2)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.text = element_text(size = rel(0.5)), legend.title = element_blank()) +
  guides(fill = guide_legend(nrow = 4)) +
  xlab("Library") +
  ylab("Number of expressed genes")

p2 <- ggplot(exp_genes, aes(x = Var2, y = raw_expr, fill = Var2)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.text = element_text(size = rel(0.5)), legend.title = element_blank()) +
  guides(fill = guide_legend(nrow = 4)) +
  xlab("Library") +
  ylab("Total raw expression")

log2abund_forplot <- reshape2::melt(as.matrix(txi$abundance)) %>%
  group_by(Var2) %>%
  mutate("log2abund" = log2(value)) %>%
  filter(!is.infinite(log2abund)) %>%
  ungroup()

p3 <- ggplot(log2abund_forplot, aes(x = Var2, y = log2abund, fill = Var2)) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "bottom", legend.text = element_text(size = rel(0.5)), legend.title = element_blank()) +
  guides(fill = guide_legend(nrow = 4)) +
  xlab("Library") +
  ylab("log2Abundance")

mu <- log2abund_forplot %>%
  group_by(Var2) %>%
  summarise(grp.mean = mean(log2abund))

p4 <- ggplot(log2abund_forplot, aes(x = log2abund, color = Var2, fill = Var2)) +
  #  geom_histogram(aes(y=..density..), position="identity", alpha=0.3) + # to much overlap, not very clear to see what's what
  geom_density(alpha = 0.3) +
  geom_vline(
    data = mu, aes(xintercept = grp.mean, color = Var2),
    linetype = "dashed", size = 1.5, alpha = 0.5
  ) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.text = element_text(size = rel(0.5)), legend.title = element_blank()) +
  guides(fill = guide_legend(nrow = 4), color = guide_legend(nrow = 4)) +
  xlab("log2Abundance")

pdf(paste(resdir, "exp-genes.pdf", sep = "/"))
  print(p)
  print(p2)
  print(p3)
  print(p4)
dev.off()

condsToCompare <- c(compcond1, compcond2)

# Decision of replicates true/false
if (repl == FALSE) { # We don't have replicates so we have to cheat - edgeR manual - 2.12 What to do if you have no replicates
  print("We don't have replicates, cheating dispersion with bcv = 0.4")
  bcv <- 0.4
  et <- exactTest(y, dispersion = bcv^2, pair = c(condsToCompare[1], condsToCompare[2]))
  d <- y # Copy for compatibility reasons downstream
  test.save <- "exactTest"
  mydesign <- paste("No replicates, used", test.save, "to test for expression differences.")
} else if (repl == TRUE) { # If we would have replicates we could do the proper testing
  print("We have replicates, hurray!")
  if ((batchCor == TRUE) && (sum(table(coldata$batch) >= 2) >= 2)) {
    print("Using batch correction.")
    mydesign <- model.matrix(~ batch + condition, data = coldata) # RK could be considered treatment of Het, we'll use intercept; otherwise change to ~0+condition
  } else {
    print("Not using batch correction.")
    mydesign <- model.matrix(~condition, data = coldata) # RK could be considered treatment of Het, we'll use intercept; otherwise change to ~0+condition
  }

  d <- estimateGLMCommonDisp(y, mydesign) # Calculate GLM for common dispersion
  d <- estimateGLMTrendedDisp(d, mydesign) # Calculate GLM for trended dispersion
  d <- estimateGLMTagwiseDisp(d, mydesign) # Calculate GLM for tagwise dispersion

  if (sum(table(coldata$condition)[names(table(coldata$condition)) %in% c(condsToCompare[1], condsToCompare[2])] > 2) >= 2) {
    print("We have OK number of replicates, going for edgeR::glmFit & edgeR::glmLRT. If you have >=6 replicates consider switching to DESeq2.")
    test.save <- "likelihood" # help for next calculations

    fit_tgw <- glmFit(d, mydesign, dispersion = d$tagwise.dispersion)
  } else {
    print("We have low number of replicates, going for edgeR::glmQLFit & edgeR::glmQLFTest.")
    test.save <- "quasi-likelihood" # help for next calculations

    fit_tgw <- glmQLFit(d, mydesign, dispersion = d$tagwise.dispersion) # fit_tgw<-glmFit(d, mydesign, dispersion=d$tagwise.dispersion); Fit tagwise dispersion;  fit_tgw<-glmQLFit(d, mydesign, dispersion=d$tagwise.dispersion) can be used if the number of replicates is low; QL (glmQLFit and glmQLFTest) is more strict in the assumptions ~ increases adj.pvalues
  }

  if (!length(grep(paste0("condition", condsToCompare[1], "$"), colnames(fit_tgw$design)))) { # If I cannot find condsToCompare[1] (usually intercept) set contrast as 1 for condsToCompare2
    my.contrasts <- makeContrasts( # Should create contrasts; if we have intercept and we want to compare coef=2 it should be the same as contrast=c(0, -1, 0) but it might be misunderstood because there is actually no contrast https://www.biostars.org/p/102036/
      postvspre = paste0(colnames(fit_tgw$design)[grep(paste0("condition", condsToCompare[2], "$"), colnames(fit_tgw$design))]), levels = fit_tgw$design
    ) # https://stackoverflow.com/questions/26813667/how-to-use-grep-to-find-exact-match
  } else { # Else make proper contrast
    my.contrasts <- makeContrasts( # Should create contrasts; if we have intercept and we want to compare coef=2 it should be the same as contrast=c(0, -1, 0) but it might be misunderstood because there is actually no contrast https://www.biostars.org/p/102036/
      postvspre = paste0(
        colnames(fit_tgw$design)[grep(paste0("condition", condsToCompare[2], "$"), colnames(fit_tgw$design))],
        "-",
        colnames(fit_tgw$design)[grep(paste0("condition", condsToCompare[1], "$"), colnames(fit_tgw$design))]
      ), levels = fit_tgw$design
    )
  }
  colnames(my.contrasts) <- "contrast"

  if (test.save == "likelihood") {
    et <- glmLRT(fit_tgw, contrast = my.contrasts[, "contrast"]) # If we have 3 conditions and want to compare 3 vs 1 we set contrast=c(0, -1, 1), if we want to compare 3 vs 1 or 2 vs 1 we just set coef=3 or coef=2, respectively; some more examples of contrast https://www.biostars.org/p/110861/
  } else if (test.save == "quasi-likelihood") {
    et <- glmQLFTest(fit_tgw, contrast = my.contrasts[, "contrast"]) # If we have 3 conditions and want to compare 3 vs 1 we set contrast=c(0, -1, 1), if we want to compare 3 vs 1 or 2 vs 1 we just set coef=3 or coef=2, respectively; some more examples of contrast https://www.biostars.org/p/110861/
  }
} else {
  stop("Don't know what to do, --replicates not set to TRUE or FALSE")
}

sink(paste(resdir, "design-control.txt", sep = "/"))
  subset(cbind(coldata, d$samples), select = -c(group))
sink()

sink(paste(resdir, "formula-control.txt", sep = "/"))
  mydesign
sink()

sink(paste(resdir, "de-summary.txt", sep = "/"))
  print(paste("Compared conditions:", paste(condsToCompare, collapse = ":")))
  print(paste("Used test:", test.save))
  print(paste("LogFC:", fcset))
  print(paste("Adj.p-value:", pvalset))
  print(summary(decideTestsDGE(et, adjust.method = "BH", p.value = pvalset, lfc = fcset))) # Summary
sink()

et.adj <- topTags(et, n = Inf, adjust.method = "BH", sort.by = "p.value") # Add FDR adjustment/correction

# Do some simple plots
pdf(paste(resdir, "de-plot.pdf", sep = "/"))
  red <- et.adj$table$logFC >= fcset & et.adj$table$FDR < pvalset
  blue <- et.adj$table$logFC <= (-fcset) & et.adj$table$FDR < pvalset
  col <- rep("black", nrow(et.adj))
  col[red] <- "red"
  col[blue] <- "blue"
  size <- rep(0.5, nrow(et.adj))
  size[red] <- 1
  size[blue] <- 1
  
  plot(et.adj$table$logCPM, et.adj$table$logFC,
    pch = 20,
    col = col, cex = size, xlab = "Average logCPM", ylab = "logFC",
    main = paste(condsToCompare[2], "vs", condsToCompare[1])
  )
  abline(h = c(-fcset, fcset), col = "blue")
  legend("topright",
    legend = c(
      paste0("NotSig", " (", sum(col == "black"), ")"),
      paste0("Up", " (", sum(col == "red"), ")"),
      paste0("Down", " (", sum(col == "blue"), ")")
    ),
    col = c("black", "red", "blue"), pch = 20, cex = 1
  )
  
  plot(et.adj$table$logFC, -log10(et.adj$table$FDR),
    pch = 20,
    col = col, cex = size, xlab = "logFC", ylab = "-log10(adj.p-val)",
    main = paste(condsToCompare[2], "vs", condsToCompare[1])
  )
  abline(v = c(-fcset, fcset))
  legend("topright",
    legend = c(
      paste0("NotSig", " (", sum(col == "black"), ")"),
      paste0("Up", " (", sum(col == "red"), ")"),
      paste0("Down", " (", sum(col == "blue"), ")")
    ),
    col = c("black", "red", "blue"), pch = 20, cex = 1
  )
dev.off()

et <- as.data.frame(et.adj)
et <- et %>%
  dplyr::rename(PValueAdj = "FDR")
et$gene_id <- rownames(et)
cts <- as.data.frame(cts)
colnames(cts) <- paste0(colnames(cts), "-counts")
et <- et %>%
  left_join(cts %>%
    mutate(gene_id = rownames(cts)))
abund.tmp <- as.data.frame(txi$abundance)
abund.tmp <- abund.tmp[, coldata$name] # Make sure abund.tmp and design are ordered the same
colnames(abund.tmp) <- paste0(colnames(abund.tmp), "-abundance")
abund.tmp$gene_id <- rownames(abund.tmp)
et <- et %>%
  left_join(abund.tmp) %>%
  left_join(desc)

# Add abundance average = average of TPM https://support.bioconductor.org/p/84883/
if (length(grep("-abundance", colnames(et))) > 0) {
  print("Using abundance from the input counts.")
} else {
  print("Using abundance from CPM calculated by edgeR")
  d2 <- d
  d2 <- calcNormFactors(d2)
  cpmabund <- cpm(d2, log = FALSE, prior.count = 1)

  cpmabund <- as.data.frame(cpmabund)
  colnames(cpmabund) <- paste0(colnames(cpmabund), "-abundance")
  cpmabund$gene_id <- rownames(cpmabund)

  et <- et %>%
    left_join(cpmabund)
}

et$logAbundance <- et %>%
  dplyr::select(ends_with("-abundance")) %>%
  rowMeans() %>%
  log2()

# Add expression bins
# Make 20 even groups by expression
quants <- quantile(et$logAbundance, probs = seq(0, 1, 0.05)) # just print the quartiles
quants
if (any(duplicated(quants))) { # If we have any duplicated quants just add a tiny number
  print("There are duplicated quantile limits, adding a tiny number.")
  quants[duplicated(quants)] <- quants[duplicated(quants)] + runif(sum(duplicated(quants)), 1e-9, 1e-8)
}
et <- within(et, AbundanceBin <- as.integer(cut(et$logAbundance, quants, include.lowest = TRUE))) # https://stackoverflow.com/questions/7508229/how-to-create-a-column-with-a-quartile-rank

et <- et[order(et$`PValueAdj`, et$`PValue`, -abs(et$`logFC`), -et$`logAbundance`), ]

# Reorganize columns
abund_cols <- grep("Abundance", colnames(et))
cols <- c(1:grep("gene_id", colnames(et)), (ncol(et) - 2):(ncol(et)) - length(abund_cols))
cols <- c(abund_cols, cols)
cols <- c(cols, (1:ncol(et))[!(1:ncol(et) %in% cols)])
et <- et[, cols]

write.table(et, paste(resdir, "de-table.xls", sep = "/"),
  quote = FALSE, sep = "\t", row.names = F
)