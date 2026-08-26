#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(VennDiagram)
  library(grid)
  library(edgeR)
  library(limma)
  library(yaml)
})

ds <- readRDS(snakemake@input[["deseq2"]])
eg <- readRDS(snakemake@input[["edger"]])
lm <- readRDS(snakemake@input[["limma"]])
diag <- readRDS(snakemake@input[["diagnostics"]])
cfg <- yaml::read_yaml(snakemake@input[["config"]])
lfc_cut <- cfg$deseq2$log2fc_threshold
species <- c("experimentalis", "gadabouti")

named_outputs <- function(key) setNames(as.character(snakemake@output[[key]]), species)
venn_out <- named_outputs("venn")
density_out <- named_outputs("density")
bcv_out <- named_outputs("bcv")
voom_out <- named_outputs("voom")
ma_paths <- as.character(snakemake@output[["ma"]])

sig_status <- function(sig, lfc) case_when(
  !is.na(sig) & sig == "yes" & lfc > 0 ~ "Up",
  !is.na(sig) & sig == "yes" & lfc < 0 ~ "Down",
  TRUE ~ "Not significant"
)
ma_plot <- function(df, abundance, lfc, sig, method, normalization, sp, outfile) {
  x <- df[[abundance]]
  if (method == "DESeq2") x <- log10(x + 1)
  pdat <- data.frame(abundance = x, logFC = df[[lfc]],
                     status = sig_status(df[[sig]], df[[lfc]]))
  p <- ggplot(pdat, aes(abundance, logFC, color = status)) +
    geom_point(alpha = .45, size = .6, na.rm = TRUE) +
    geom_hline(yintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    scale_color_manual(values = c(Up = "#E64B35", Down = "#4DBBD5",
                                  `Not significant` = "grey75")) +
    labs(title = paste0("MA plot — ", method, " (", sp, ")"),
         subtitle = paste("Normalization:", normalization),
         x = if (method == "DESeq2") "log10(mean normalized count + 1)" else abundance,
         y = "log2 fold change", color = NULL) + theme_classic(base_size = 11)
  ggsave(outfile, p, width = 7, height = 5)
}

for (sp in species) {
  sets <- list(DESeq2 = ds[[sp]]$GeneID[ds[[sp]]$sig == "yes" & !is.na(ds[[sp]]$sig)],
               edgeR = eg[[sp]]$GeneID[eg[[sp]]$sig == "yes" & !is.na(eg[[sp]]$sig)],
               limma = lm[[sp]]$GeneID[lm[[sp]]$sig == "yes" & !is.na(lm[[sp]]$sig)])
  vp <- venn.diagram(sets, filename = NULL,
    category.names = c("DESeq2", "edgeR", "limma"),
    main = paste("Overlap of significant genes:", sp), cex = 1.3,
    cat.cex = 1.2, main.cex = 1.5)
  pdf(venn_out[[sp]], width = 7, height = 7); grid.draw(vp); dev.off()

  den <- bind_rows(
    transmute(ds[[sp]], Log2FC = log2FoldChange, Method = "DESeq2"),
    transmute(eg[[sp]], Log2FC = logFC, Method = "edgeR"),
    transmute(lm[[sp]], Log2FC = logFC, Method = "limma")
  ) %>% filter(is.finite(Log2FC))
  p <- ggplot(den, aes(Log2FC, fill = Method)) + geom_density(alpha = .5) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    labs(title = paste("Log2FC Distribution:", sp), x = "log2(Fold Change)", y = "Density") +
    theme_minimal()
  ggsave(density_out[[sp]], p, width = 7, height = 5)

  pdf(bcv_out[[sp]], width = 7, height = 5)
  plotBCV(diag[[sp]]$dge, main = paste("edgeR Dispersion:", sp))
  dev.off()

  pdf(voom_out[[sp]], width = 7, height = 5)
  voom(diag[[sp]]$dge, diag[[sp]]$design, plot = TRUE)
  title(main = paste("voom Mean-Variance Trend:", sp))
  dev.off()

  ma_plot(ds[[sp]], "baseMean", "log2FoldChange", "sig", "DESeq2",
          "median-of-ratios", sp,
          ma_paths[grepl(paste0("MA_DESeq2_", sp, "\\.pdf$"), ma_paths)])
  ma_plot(eg[[sp]], "logCPM", "logFC", "sig", "edgeR", "TMM", sp,
          ma_paths[grepl(paste0("MA_edgeR_", sp, "\\.pdf$"), ma_paths)])
  ma_plot(lm[[sp]], "AveExpr", "logFC", "sig", "limma-voom", "TMM + voom", sp,
          ma_paths[grepl(paste0("MA_limma_voom_", sp, "\\.pdf$"), ma_paths)])
}

