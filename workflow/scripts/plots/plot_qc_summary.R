#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

qc <- read.csv(snakemake@input[["qc_summary"]])
colors_condition <- c(control = "#1b9e77", anhydrobiosis = "#d95f02")
set.seed(123)

p1 <- ggplot(qc, aes(species, million_reads, fill = condition)) +
  geom_boxplot() + geom_jitter(width = 0.2, size = 2) +
  scale_fill_manual(values = colors_condition) +
  labs(title = "Million Reads per Sample") + theme_minimal(base_size = 12)
p2 <- ggplot(qc, aes(species, features_above_threshold, fill = condition)) +
  geom_boxplot() + geom_jitter(width = 0.2, size = 2) +
  scale_fill_manual(values = colors_condition) +
  labs(title = "Features Above Threshold") + theme_minimal(base_size = 12)
p3 <- ggplot(qc, aes(species, percent_aligned, fill = condition)) +
  geom_boxplot() + geom_jitter(width = 0.2, size = 2) +
  scale_fill_manual(values = colors_condition) +
  labs(title = "% Aligned") + theme_minimal(base_size = 12)

ggsave(snakemake@output[["plot"]], p3 + p1 + p2 + plot_layout(ncol = 3),
       width = 18, height = 5, dpi = 300)

