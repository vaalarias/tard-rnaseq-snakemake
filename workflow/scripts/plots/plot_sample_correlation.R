#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
})

colors_condition <- c(control = "#1b9e77", anhydrobiosis = "#d95f02")
colors_species <- c(experimentalis = "#56b4e9", gadabouti = "#4b2588ff")
vst <- readRDS(snakemake@input[["vst"]])
raw <- readRDS(snakemake@input[["raw_counts"]])
mat <- if (is.matrix(vst)) vst else assay(vst)
coldata <- raw$coldata
rownames(coldata) <- coldata$sample

cor_matrix <- cor(mat)
cor_no_diag <- cor_matrix
diag(cor_no_diag) <- NA
min_cor <- min(cor_no_diag, na.rm = TRUE)
max_cor <- max(cor_no_diag, na.rm = TRUE)
mid_point <- (min_cor + max_cor) / 2
col_fun <- colorRamp2(c(min_cor, mid_point, 1), c("#313695", "#ffffff", "#a50026"))

ann <- coldata[colnames(cor_matrix), c("condition", "species"), drop = FALSE]
ht <- Heatmap(
  cor_matrix, name = "Correlation", col = col_fun,
  top_annotation = HeatmapAnnotation(
    df = ann,
    col = list(condition = colors_condition, species = colors_species)
  ),
  show_row_names = TRUE, show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 9), column_names_rot = 45,
  cluster_rows = TRUE, cluster_columns = TRUE,
  heatmap_legend_param = list(at = c(.6, .8, 1), labels = c(.6, .8, 1))
)
pdf(snakemake@output[["plot"]], width = 8, height = 8)
draw(ht)
dev.off()

