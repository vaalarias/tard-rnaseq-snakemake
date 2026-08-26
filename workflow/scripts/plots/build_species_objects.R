#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(DESeq2))

raw <- readRDS(snakemake@input[["raw_counts"]])
counts <- raw$counts
coldata <- raw$coldata

if (!"sample" %in% colnames(coldata)) stop("coldata requires a 'sample' column")
rownames(coldata) <- coldata$sample
if (!identical(colnames(counts), rownames(coldata))) {
  coldata <- coldata[colnames(counts), , drop = FALSE]
}

dds_list <- list()
vst_list <- list()
for (sp in unique(coldata$species)) {
  samples_sp <- rownames(coldata)[coldata$species == sp]
  dds <- DESeqDataSetFromMatrix(
    countData = counts[, samples_sp, drop = FALSE],
    colData = coldata[samples_sp, , drop = FALSE],
    design = ~ condition
  )
  dds <- DESeq(dds)
  dds_list[[sp]] <- dds
  vst_list[[sp]] <- vst(dds)
}

saveRDS(dds_list, snakemake@output[["dds_list"]])
saveRDS(vst_list, snakemake@output[["vst_list"]])

