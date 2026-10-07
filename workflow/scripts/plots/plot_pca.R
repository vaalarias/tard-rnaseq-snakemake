#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(matrixStats)
    library(patchwork)
    library(yaml)
})

colors_condition <- c(
    control = "#1b9e77",
    anhydrobiosis = "#d95f02"
)

shapes_species <- c(
    experimentalis = 16,
    gadabouti = 17
)

condition_labels <- c(
    control = "Control",
    anhydrobiosis = "Post-anhydrobiosis"
)

species_labels <- c(
    experimentalis = "Pam. experimentalis",
    gadabouti = "Pam. gadabouti"
)

cm_to_in <- function(x) x / 2.54

cfg <- yaml::read_yaml(
    snakemake@input[["plot_metadata"]]
)

raw <- readRDS(
    snakemake@input[["raw_counts"]]
)

coldata <- raw$coldata
rownames(coldata) <- coldata$sample

to_matrix <- function(x) {
    if (is.matrix(x)) {
        x
    } else {
        assay(x)
    }
}

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

pca_theme <- theme_minimal(base_size = 10) +
    theme(
        plot.title = element_text(
            face = "bold",
            hjust = 0.5,
            size = 11
        ),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.title = element_text(
            face = "bold",
            size = 9
        ),
        legend.text = element_text(size = 8),
        plot.margin = margin(6, 8, 6, 8)
    )

# Combined PCA
vst_all <- readRDS(
    snakemake@input[["vst"]]
)

mat_all <- to_matrix(vst_all)

pc_all <- make_pca_data(
    mat_all,
    coldata[colnames(mat_all), , drop = FALSE],
    cfg$pca_ntop
)

p_all <- ggplot(
    pc_all$data,
    aes(
        x = PC1,
        y = PC2,
        color = condition,
        shape = species
    )
) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(
        values = colors_condition,
        breaks = names(condition_labels),
        labels = condition_labels
    ) +
    scale_shape_manual(
        values = shapes_species,
        breaks = names(species_labels),
        labels = species_labels
    ) +
    labs(
        title = "All libraries",
        x = paste0("PC1 (", pc_all$pct[1], "%)"),
        y = paste0("PC2 (", pc_all$pct[2], "%)"),
        color = "Condition",
        shape = "Species"
    ) +
    pca_theme

# Species-specific PCA
vst_list <- readRDS(
    snakemake@input[["vst_list"]]
)

make_species_plot <- function(sp, title) {

    mat <- to_matrix(vst_list[[sp]])

    meta <- coldata[
        colnames(mat),
        ,
        drop = FALSE
    ]

    pc <- make_pca_data(
        mat,
        meta,
        cfg$pca_ntop
    )

    ggplot(
        pc$data,
        aes(
            x = PC1,
            y = PC2,
            color = condition
        )
    ) +
        geom_point(size = 3, alpha = 0.9) +
        scale_color_manual(
            values = colors_condition,
            breaks = names(condition_labels),
            labels = condition_labels
        ) +
        labs(
            title = title,
            x = paste0("PC1 (", pc$pct[1], "%)"),
            y = paste0("PC2 (", pc$pct[2], "%)"),
            color = "Condition"
        ) +
        pca_theme
}

p_experimentalis <- make_species_plot(
    "experimentalis",
    "Pam. experimentalis"
)

p_gadabouti <- make_species_plot(
    "gadabouti",
    "Pam. gadabouti"
)

# Top: combined PCA
# Bottom: species-specific PCA
pca_panel <- (
    p_all /
        (p_experimentalis | p_gadabouti) +
        plot_layout(
            heights = c(1.15, 1),
            guides = "collect"
        ) +
        plot_annotation(
            tag_levels = "A"
        )
) &
    theme(
        legend.position = "bottom",
        plot.tag = element_text(
            face = "bold",
            size = 13
        )
    )

ggsave(
    filename = snakemake@output[["plot"]],
    plot = pca_panel,
    width = cm_to_in(18),
    height = cm_to_in(15),
    dpi = 300,
    device = cairo_pdf
)
