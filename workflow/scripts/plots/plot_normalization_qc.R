#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

raw <- readRDS(snakemake@input[["raw_counts"]])
tpm <- as.matrix(readRDS(snakemake@input[["tpm"]]))
vst_obj <- readRDS(snakemake@input[["vst"]])
vst <- if (is.matrix(vst_obj)) vst_obj else assay(vst_obj)
coldata <- raw$coldata
colors_species <- c(experimentalis = "#56b4e9", gadabouti = "#4b2588ff")

tpm_long <- as.data.frame(tpm) %>% rownames_to_column("gene_id") %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "TPM") %>%
  left_join(coldata[, c("sample", "species", "condition")], by = "sample") %>%
  mutate(logTPM = log10(TPM + 1))
vst_long <- as.data.frame(vst) %>% rownames_to_column("gene_id") %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "VST") %>%
  left_join(coldata[, c("sample", "species", "condition")], by = "sample")

p1 <- ggplot(tpm_long, aes(logTPM)) +
  geom_histogram(bins = 100, color = "black", fill = "#4C78A8") +
  facet_wrap(~species, scales = "free_y") + theme_bw() +
  labs(title = "TPM distribution", x = expression(log[10](TPM+1)), y = "Gene count")
p2 <- ggplot(tpm_long, aes(sample, logTPM, fill = species)) +
  geom_boxplot(outlier.size = .15) + scale_fill_manual(values = colors_species) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "TPM distribution per sample", x = "", y = expression(log[10](TPM+1)))
p3 <- ggplot(vst_long, aes(VST)) +
  geom_histogram(bins = 100, color = "black", fill = "#F58518") +
  facet_wrap(~species, scales = "free_y") + theme_bw() +
  labs(title = "VST distribution", x = "Variance Stabilizing Transformation", y = "Gene count")
p4 <- ggplot(vst_long, aes(sample, VST, fill = species)) +
  geom_boxplot(outlier.size = .15) + scale_fill_manual(values = colors_species) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "VST distribution per sample", x = "", y = "VST")

ggsave(snakemake@output[["tpm_hist"]], p1, width = 8, height = 5)
ggsave(snakemake@output[["tpm_box"]], p2, width = 9, height = 5)
ggsave(snakemake@output[["vst_hist"]], p3, width = 8, height = 5)
ggsave(snakemake@output[["vst_box"]], p4, width = 9, height = 5)

