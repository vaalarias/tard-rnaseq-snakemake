#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

res_list <- readRDS(snakemake@input[["deseq_results"]])
outputs <- setNames(as.character(snakemake@output), c("experimentalis", "gadabouti"))
display_names <- c(experimentalis = "Experimentalis", gadabouti = "Gadabouti")
sig_colors <- c(Up = "#E64B35", Down = "#4DBBD5", `Not Sig` = "grey70")

for (sp in names(outputs)) {
  res <- res_list[[sp]] %>% mutate(sig_status = case_when(
    sig == "yes" & log2FoldChange > 0 ~ "Up",
    sig == "yes" & log2FoldChange < 0 ~ "Down",
    TRUE ~ "Not Sig"
  ))
  n_up <- sum(res$sig_status == "Up", na.rm = TRUE)
  n_down <- sum(res$sig_status == "Down", na.rm = TRUE)
  finite_y <- -log10(res$padj[is.finite(-log10(res$padj))])
  x_min <- min(res$log2FoldChange, na.rm = TRUE)
  x_max <- max(res$log2FoldChange, na.rm = TRUE)
  y_min <- min(finite_y, na.rm = TRUE)
  y_max <- max(finite_y, na.rm = TRUE)
  p <- ggplot(res, aes(log2FoldChange, -log10(padj), color = sig_status)) +
    geom_point(alpha = 0.6, size = 2) + scale_color_manual(values = sig_colors) +
    geom_vline(xintercept = c(-1.6, 1.6), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    labs(title = paste0("Volcano Plot - ", display_names[[sp]]),
         x = "log2 Fold Change", y = "-log10(adj. p-value)", color = "Significance") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "right", plot.title = element_text(face = "bold", size = 16)) +
    annotate("text", x = x_max, y = y_max * .95, label = paste0("Up: ", n_up),
             hjust = 1, vjust = 1, size = 5, color = "grey40") +
    annotate("text", x = x_min, y = y_max * .95, label = paste0("Down: ", n_down),
             hjust = 0, vjust = 1, size = 5, color = "grey40") +
    annotate("text", x = x_max, y = y_min, label = paste0("Total sig: ", n_up + n_down),
             hjust = 1, vjust = 0, size = 5, color = "black")
  ggsave(outputs[[sp]], p, width = 7, height = 6)
}

