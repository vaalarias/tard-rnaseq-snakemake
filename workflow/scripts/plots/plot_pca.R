#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(matrixStats)
  library(yaml)
})

colors_condition <- c(control = "#1b9e77", anhydrobiosis = "#d95f02")
colors_species <- c(experimentalis = "#56b4e9", gadabouti = "#4b2588ff")
cm_to_in <- function(x) x / 2.54
cfg <- yaml::read_yaml(snakemake@input[["plot_metadata"]])
raw <- readRDS(snakemake@input[["raw_counts"]])
coldata <- raw$coldata
rownames(coldata) <- coldata$sample

to_matrix <- function(x) if (is.matrix(x)) x else assay(x)
make_pca_data <- function(mat, meta, ntop = 500) {
  mat <- as.matrix(mat)

  rv <- matrixStats::rowVars(
    mat,
    useNames = TRUE
  )

  keep <- is.finite(rv) & rv > 0

  mat <- mat[keep, , drop = FALSE]
  rv <- rv[keep]

  if (nrow(mat) < 2) {
    stop(
      "PCA requires at least two genes with non-zero variance."
    )
  }

  ntop <- min(
    as.integer(ntop),
    length(rv)
  )

  selected <- order(
    rv,
    decreasing = TRUE
  )[seq_len(ntop)]

  pc <- prcomp(
    t(mat[selected, , drop = FALSE]),
    center = TRUE,
    scale. = FALSE
  )

  pct <- round(
    100 * pc$sdev^2 / sum(pc$sdev^2),
    1
  )

  metadata <- meta[
    rownames(pc$x),
    ,
    drop = FALSE
  ]

  df <- cbind(
    as.data.frame(
      pc$x[, 1:2, drop = FALSE]
    ),
    metadata
  )

  list(
    data = df,
    pct = pct
  )
}

# Total PCA, equivalent to DESeq2::plotPCA(ntop=500).
vst_all <- readRDS(snakemake@input[["vst"]])
mat_all <- to_matrix(vst_all)
pc <- make_pca_data(mat_all, coldata[colnames(mat_all), , drop = FALSE], cfg$pca_ntop)
p <- ggplot(pc$data, aes(PC1, PC2, color = condition, shape = species)) +
  geom_point(size = 2) + scale_color_manual(values = colors_condition) +
  labs(title = "PCA (VST Normalization)",
       x = paste0("PC1: ", pc$pct[1], "%"), y = paste0("PC2: ", pc$pct[2], "%")) +
  theme_minimal(base_size = 6)
ggsave(snakemake@output[["total"]], p, width = cm_to_in(8), height = cm_to_in(5), dpi = 300)

# Per-species PCA from the species-specific VST objects, as in the notebook.
vst_list <- readRDS(snakemake@input[["vst_list"]])
out_by_species <- setNames(as.character(snakemake@output[["species"]]),
                           c("experimentalis", "gadabouti"))
for (sp in names(out_by_species)) {
  mat <- to_matrix(vst_list[[sp]])
  meta <- coldata[colnames(mat), , drop = FALSE]
  pc <- make_pca_data(mat, meta, cfg$pca_ntop)
  p <- ggplot(pc$data, aes(PC1, PC2, color = condition)) +
    geom_point(size = 2) + scale_color_manual(values = colors_condition) +
    ggtitle(paste("PCA -", sp)) +
    xlab(paste0("PC1 (", round(pc$pct[1]), "%)")) +
    ylab(paste0("PC2 (", round(pc$pct[2]), "%)")) +
    theme_minimal(base_size = 6)
  ggsave(out_by_species[[sp]], p, width = cm_to_in(8), height = cm_to_in(5), dpi = 300)
}

