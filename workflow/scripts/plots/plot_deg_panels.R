#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(dplyr)
    library(ggplot2)
    library(cowplot)
    library(grid)
})

# ============================================================
# Colours and labels
# ============================================================

colors_condition <- c(
    control = "#1b9e77",
    anhydrobiosis = "#d95f02"
)

condition_labels <- c(
    control = "Control",
    anhydrobiosis = "Post-anhydrobiosis"
)

sig_colors <- c(
    Up = "#E64B35",
    Down = "#4DBBD5",
    `Not Sig` = "grey70"
)

display_names <- c(
    experimentalis = "Experimentalis",
    gadabouti = "Gadabouti"
)

required_species <- c(
    "experimentalis",
    "gadabouti"
)

# ============================================================
# Helper functions
# ============================================================

to_matrix <- function(x) {
    if (is.matrix(x)) {
        return(x)
    }

    as.matrix(
        SummarizedExperiment::assay(x)
    )
}

make_volcano <- function(results, species_name) {

    required_columns <- c(
        "GeneID",
        "log2FoldChange",
        "padj",
        "sig"
    )

    missing_columns <- setdiff(
        required_columns,
        colnames(results)
    )

    if (length(missing_columns)) {
        stop(
            "Missing columns in DESeq2 results for ",
            species_name,
            ": ",
            paste(missing_columns, collapse = ", ")
        )
    }

    res <- results %>%
        mutate(
            sig_status = case_when(
                sig == "yes" &
                    log2FoldChange > 0 ~ "Up",

                sig == "yes" &
                    log2FoldChange < 0 ~ "Down",

                TRUE ~ "Not Sig"
            ),
            sig_status = factor(
                sig_status,
                levels = c(
                    "Up",
                    "Down",
                    "Not Sig"
                )
            )
        )

    n_up <- sum(
        res$sig_status == "Up",
        na.rm = TRUE
    )

    n_down <- sum(
        res$sig_status == "Down",
        na.rm = TRUE
    )

    finite_lfc <- res$log2FoldChange[
        is.finite(res$log2FoldChange)
    ]

    transformed_padj <- -log10(
        res$padj
    )

    finite_y <- transformed_padj[
        is.finite(transformed_padj)
    ]

    if (!length(finite_lfc)) {
        stop(
            "No finite log2 fold-change values found for ",
            species_name
        )
    }

    if (!length(finite_y)) {
        stop(
            "No finite adjusted p-values found for ",
            species_name
        )
    }

    x_min <- min(finite_lfc)
    x_max <- max(finite_lfc)
    y_min <- min(finite_y)
    y_max <- max(finite_y)

    ggplot(
        res,
        aes(
            x = log2FoldChange,
            y = -log10(padj),
            color = sig_status
        )
    ) +
        geom_point(
            alpha = 0.6,
            size = 2,
            na.rm = TRUE
        ) +
        scale_color_manual(
            values = sig_colors,
            drop = FALSE
        ) +
        geom_vline(
            xintercept = c(-1.6, 1.6),
            linetype = "dashed",
            color = "black"
        ) +
        geom_hline(
            yintercept = -log10(0.05),
            linetype = "dashed",
            color = "black"
        ) +
        labs(
            title = paste0(
                "Pam. ",
                display_names[[species_name]]
            ),
            x = "log2 Fold Change",
            y = "-log10(adj. p-value)",
            color = "Significance"
        ) +
        theme_minimal(
            base_size = 14
        ) +
        theme(
            legend.position = "right",
            plot.title = element_text(
                face = "bold",
                size = 16
            ),
            panel.grid.minor = element_blank()
        ) +
        annotate(
            "text",
            x = x_max,
            y = y_max * 0.95,
            label = paste0(
                "Up: ",
                n_up
            ),
            hjust = 1,
            vjust = 1,
            size = 5,
            color = "grey40"
        ) +
        annotate(
            "text",
            x = x_min,
            y = y_max * 0.95,
            label = paste0(
                "Down: ",
                n_down
            ),
            hjust = 0,
            vjust = 1,
            size = 5,
            color = "grey40"
        ) 
}
make_heatmap <- function(
    z,
    annotation_data,
    heatmap_col_fun,
    species_name
) {
    Heatmap(
        z,
        name = "Z-score",
        col = heatmap_col_fun,

        top_annotation = HeatmapAnnotation(
            df = annotation_data,
            col = list(
                condition = colors_condition
            ),
            annotation_legend_param = list(
                condition = list(
                    title = "Condition",
                    at = names(condition_labels),
                    labels = unname(condition_labels)
                )
            )
        ),

        show_row_names = FALSE,
        show_column_names = FALSE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,

        column_title = paste(
            "Pam.",
            display_names[[species_name]]
        ),
        column_title_gp = gpar(
            fontsize = 16,
            fontface = "bold"
        ),

        column_names_gp = gpar(
            fontsize = 9
        ),

        heatmap_legend_param = list(
            title = "Z-score",
            labels_gp = gpar(fontsize = 8),
            title_gp = gpar(
                fontsize = 9,
                fontface = "bold"
            )
        )
    )
}



# ============================================================
# Read inputs
# ============================================================

res_list <- readRDS(
    snakemake@input[["deseq_results"]]
)

vst_list <- readRDS(
    snakemake@input[["vst_list"]]
)

raw <- readRDS(
    snakemake@input[["raw_counts"]]
)

coldata <- raw$coldata

if (!"sample" %in% colnames(coldata)) {
    stop(
        "coldata must contain a 'sample' column."
    )
}

missing_results <- setdiff(
    required_species,
    names(res_list)
)

missing_vst <- setdiff(
    required_species,
    names(vst_list)
)

if (length(missing_results)) {
    stop(
        "Species missing from DESeq2 results: ",
        paste(missing_results, collapse = ", ")
    )
}

if (length(missing_vst)) {
    stop(
        "Species missing from VST objects: ",
        paste(missing_vst, collapse = ", ")
    )
}

# ============================================================
# Prepare heatmap matrices
# ============================================================

centered <- list()
all_z_values <- numeric()

for (sp in required_species) {

    mat <- to_matrix(
        vst_list[[sp]]
    )

    species_results <- res_list[[sp]]

    sig_genes <- species_results %>%
        filter(
            !is.na(GeneID),
            !is.na(sig),
            sig == "yes"
        ) %>%
        pull(GeneID) %>%
        unique()

    sig_genes <- intersect(
        sig_genes,
        rownames(mat)
    )

    if (!length(sig_genes)) {
        stop(
            "No significant genes overlap the VST matrix for ",
            sp
        )
    }

    expression_matrix <- mat[
        sig_genes,
        ,
        drop = FALSE
    ]

    gene_sd <- apply(
        expression_matrix,
        1,
        sd,
        na.rm = TRUE
    )

    keep <- is.finite(gene_sd) &
        gene_sd > 0

    expression_matrix <- expression_matrix[
        keep,
        ,
        drop = FALSE
    ]

    if (!nrow(expression_matrix)) {
        stop(
            "All significant genes have zero variance for ",
            sp
        )
    }

    z <- t(
        scale(
            t(expression_matrix)
        )
    )

    z[!is.finite(z)] <- 0

    centered[[sp]] <- z

    all_z_values <- c(
        all_z_values,
        as.numeric(z)
    )

    message(
        sp,
        ": ",
        nrow(z),
        " significant genes included in the heatmap"
    )
}

# Use the same symmetric Z-score scale for both species

max_abs <- max(
    abs(all_z_values),
    na.rm = TRUE
)

if (!is.finite(max_abs) || max_abs == 0) {
    max_abs <- 1
}

heatmap_col_fun <- colorRamp2(
    c(
        -max_abs,
        0,
        max_abs
    ),
    c(
        "#313695",
        "white",
        "#a50026"
    )
)

# ============================================================
# Define named outputs
# ============================================================

volcano_outputs <- c(
    experimentalis = as.character(
        snakemake@output[[
            "volcano_experimentalis"
        ]]
    ),
    gadabouti = as.character(
        snakemake@output[[
            "volcano_gadabouti"
        ]]
    )
)

heatmap_outputs <- c(
    experimentalis = as.character(
        snakemake@output[[
            "heatmap_experimentalis"
        ]]
    ),
    gadabouti = as.character(
        snakemake@output[[
            "heatmap_gadabouti"
        ]]
    )
)

panel_outputs <- c(
    experimentalis = as.character(
        snakemake@output[[
            "panel_experimentalis"
        ]]
    ),
    gadabouti = as.character(
        snakemake@output[[
            "panel_gadabouti"
        ]]
    )
)

all_output_paths <- c(
    volcano_outputs,
    heatmap_outputs,
    panel_outputs
)

invisible(
    lapply(
        unique(dirname(all_output_paths)),
        dir.create,
        recursive = TRUE,
        showWarnings = FALSE
    )
)

# ============================================================
# Generate individual plots and panels
# ============================================================

for (sp in required_species) {

    message(
        "Generating DEG plots for ",
        sp
    )

    # --------------------------------------------------------
    # Volcano plot
    # --------------------------------------------------------

    volcano <- make_volcano(
        res_list[[sp]],
        sp
    )

    ggsave(
        filename = volcano_outputs[[sp]],
        plot = volcano,
        width = 7,
        height = 6,
        device = cairo_pdf
    )

    # --------------------------------------------------------
    # Heatmap annotation
    # --------------------------------------------------------

    z <- centered[[sp]]

    sample_match <- match(
        colnames(z),
        coldata$sample
    )

    if (anyNA(sample_match)) {
        stop(
            "Samples from the VST matrix were not found in coldata for ",
            sp,
            ": ",
            paste(
                colnames(z)[is.na(sample_match)],
                collapse = ", "
            )
        )
    }

    annotation_data <- data.frame(
        condition = as.character(
            coldata$condition[sample_match]
        ),
        row.names = colnames(z)
    )

    unknown_conditions <- setdiff(
        unique(annotation_data$condition),
        names(colors_condition)
    )

    if (length(unknown_conditions)) {
        stop(
            "Conditions without an assigned color: ",
            paste(unknown_conditions, collapse = ", ")
        )
    }

    # --------------------------------------------------------
    # Individual heatmap
    # --------------------------------------------------------

    individual_heatmap <- make_heatmap(
        z,
        annotation_data,
        heatmap_col_fun,
        sp
    )

    pdf(
        heatmap_outputs[[sp]],
        width = 8,
        height = 8,
        useDingbats = FALSE
    )

    draw(
        individual_heatmap,
        heatmap_legend_side = "right",
        annotation_legend_side = "right",
        merge_legends = TRUE
    )

    dev.off()

    # --------------------------------------------------------
    # Capture heatmap for the combined panel
    # --------------------------------------------------------

    panel_heatmap <- make_heatmap(
        z,
        annotation_data,
        heatmap_col_fun,
        sp
    )

    temporary_pdf <- tempfile(
    pattern = paste0("heatmap_", sp, "_"),
    fileext = ".pdf"
    )

    pdf(
        temporary_pdf,
        width = 8,
        height = 7,
        useDingbats = FALSE
    )

    draw(
        panel_heatmap,
        heatmap_legend_side = "right",
        annotation_legend_side = "right",
        merge_legends = TRUE
    )

    heatmap_grob <- grid.grab()

    dev.off()

    unlink(temporary_pdf)

    # --------------------------------------------------------
    # Side-by-side panel
    # --------------------------------------------------------

    panel <- plot_grid(
        volcano,
        heatmap_grob,
        ncol = 2,
        labels = c(
            "A",
            "B"
        ),
        label_fontface = "bold",
        label_size = 14,
        rel_widths = c(
            1,
            1.15
        ),
        align = "h",
        axis = "tb"
    )

    ggsave(
        filename = panel_outputs[[sp]],
        plot = panel,
        width = 15,
        height = 7,
        device = cairo_pdf
    )

    message(
        "Finished DEG plots for ",
        sp
    )
}