#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(DESeq2)
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
})

colors_condition <- c(
    control = "#1b9e77",
    anhydrobiosis = "#d95f02"
)

colors_species <- c(
    experimentalis = "#56b4e9",
    gadabouti = "#4b2588ff"
)

vst <- readRDS(
    snakemake@input[["vst"]]
)

raw <- readRDS(
    snakemake@input[["raw_counts"]]
)

mat <- if (is.matrix(vst)) {
    vst
} else {
    assay(vst)
}

coldata <- raw$coldata
rownames(coldata) <- coldata$sample

cor_matrix <- cor(
    mat,
    method = "pearson",
    use = "pairwise.complete.obs"
)

cor_no_diag <- cor_matrix
diag(cor_no_diag) <- NA

min_cor <- min(
    cor_no_diag,
    na.rm = TRUE
)

# Ensure that the displayed scale includes at least 0.6
lower_limit <- min(0.6, min_cor)
mid_point <- (lower_limit + 1) / 2

col_fun <- colorRamp2(
    c(lower_limit, mid_point, 1),
    c("#c4c6ee", "#8589c7", "#313695")
)

ann <- coldata[
    colnames(cor_matrix),
    c("condition", "species"),
    drop = FALSE
]

top_ann <- HeatmapAnnotation(
    df = ann,
    col = list(
        condition = colors_condition,
        species = colors_species
    ),
    annotation_name_gp = gpar(
        fontsize = 9,
        fontface = "bold"
    ),
    simple_anno_size = unit(0.35, "cm")
)

ht <- Heatmap(
    cor_matrix,
    name = "Correlation",
    col = col_fun,

    top_annotation = HeatmapAnnotation(
        df = ann,
        col = list(
            condition = colors_condition,
            species = colors_species
        ),
        simple_anno_size = unit(0.35, "cm")
    ),

    # Preserve the original ComplexHeatmap clustering
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows = "complete",
    clustering_method_columns = "complete",

    # Square heatmap body
    width = unit(12, "cm"),
    height = unit(12, "cm"),

    # Show sample names only on the right
    show_row_names = FALSE,
    show_column_names = FALSE,

    row_names_gp = gpar(fontsize = 8),

    border = TRUE,

    heatmap_legend_param = list(
        at = c(0.6, 0.8, 1),
        labels = c("0.6", "0.8", "1.0"),
        direction = "horizontal",
        legend_width = unit(4, "cm")
    )
)

pdf(
    snakemake@output[["plot"]],
    width = 18 / 2.54,
    height = 18 / 2.54,
    useDingbats = FALSE
)

draw(
    ht,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    merge_legends = TRUE,
    padding = unit(c(4, 4, 4, 4), "mm")
)

dev.off()

