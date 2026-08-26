#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(VennDiagram)
  library(dplyr)
  library(grid)
})

res <- readRDS(snakemake@input[["deseq_results"]])
sig_list <- list(
  Gadabouti = res[["gadabouti"]] %>% filter(!is.na(sig), sig == "yes") %>% pull(GeneID) %>% unique(),
  Experimentalis = res[["experimentalis"]] %>% filter(!is.na(sig), sig == "yes") %>% pull(GeneID) %>% unique()
)
venn_plot <- venn.diagram(
  x = sig_list, filename = NULL,
  fill = c(Gadabouti = "#4b2588ff", Experimentalis = "#56b4e9"),
  alpha = 0.5, cex = 2, fontface = "bold",
  cat.cex = 1.5, cat.fontface = "bold", cat.pos = c(-20, 20), margin = 0.1
)
pdf(snakemake@output[["plot"]], width = 6, height = 6)
grid.draw(venn_plot)
dev.off()

