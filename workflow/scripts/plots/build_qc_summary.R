#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(yaml))

raw <- readRDS(snakemake@input[["raw_counts"]])
cfg <- yaml::read_yaml(snakemake@input[["plot_metadata"]])
counts <- raw$counts
coldata <- raw$coldata
if (is.null(rownames(coldata)) || !all(colnames(counts) %in% rownames(coldata))) {
  rownames(coldata) <- coldata$sample
}

total_reads <- colSums(counts)
threshold <- cfg$features_above_count
features_above_threshold <- colSums(counts > threshold)
percent_aligned <- unlist(cfg$percent_aligned)

qc <- data.frame(
  sample = colnames(counts),
  species = coldata[colnames(counts), "species"],
  condition = coldata[colnames(counts), "condition"],
  total_reads = total_reads,
  million_reads = total_reads / 1e6,
  features_above_threshold = features_above_threshold,
  percent_aligned = as.numeric(percent_aligned[colnames(counts)]),
  row.names = NULL,
  check.names = FALSE
)

write.csv(qc, snakemake@output[["qc_summary"]], row.names = FALSE)
