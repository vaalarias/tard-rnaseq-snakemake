#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(dplyr)
    library(ggplot2)
    library(yaml)
})

# ============================================================
# Inputs
# ============================================================

annot_exp <- readRDS(
    snakemake@input[["annot_exp"]]
)

annot_gad <- readRDS(
    snakemake@input[["annot_gad"]]
)

term2gene <- readRDS(
    snakemake@input[["term2gene"]]
)

term2name <- readRDS(
    snakemake@input[["term2name"]]
)

config <- yaml::read_yaml(
    snakemake@input[["config"]]
)

# ============================================================
# Outputs
# ============================================================

pdf_outputs <- c(
    experimentalis_up = as.character(
        snakemake@output[["exp_up_pdf"]]
    ),
    experimentalis_down = as.character(
        snakemake@output[["exp_down_pdf"]]
    ),
    gadabouti_up = as.character(
        snakemake@output[["gad_up_pdf"]]
    ),
    gadabouti_down = as.character(
        snakemake@output[["gad_down_pdf"]]
    )
)

csv_outputs <- c(
    experimentalis_up = as.character(
        snakemake@output[["exp_up_csv"]]
    ),
    experimentalis_down = as.character(
        snakemake@output[["exp_down_csv"]]
    ),
    gadabouti_up = as.character(
        snakemake@output[["gad_up_csv"]]
    ),
    gadabouti_down = as.character(
        snakemake@output[["gad_down_csv"]]
    )
)

rds_outputs <- c(
    experimentalis_up = as.character(
        snakemake@output[["exp_up_rds"]]
    ),
    experimentalis_down = as.character(
        snakemake@output[["exp_down_rds"]]
    ),
    gadabouti_up = as.character(
        snakemake@output[["gad_up_rds"]]
    ),
    gadabouti_down = as.character(
        snakemake@output[["gad_down_rds"]]
    )
)

all_outputs <- c(
    pdf_outputs,
    csv_outputs,
    rds_outputs
)

invisible(
    lapply(
        unique(dirname(all_outputs)),
        dir.create,
        recursive = TRUE,
        showWarnings = FALSE
    )
)

# ============================================================
# Configuration
# ============================================================

padj_threshold <- as.numeric(
    config$deseq2$padj_threshold
)

log2fc_threshold <- as.numeric(
    config$deseq2$log2fc_threshold
)

show_category <- as.integer(
    config$plots$general$showCategory
)

plot_width <- as.numeric(
    config$plots$general$width
)

plot_height <- as.numeric(
    config$plots$general$height
)

species_colors <- unlist(
    config$plots$colors$species
)

required_species_colors <- c(
    "experimentalis",
    "gadabouti"
)

missing_species_colors <- setdiff(
    required_species_colors,
    names(species_colors)
)

if (length(missing_species_colors)) {
    stop(
        "Missing species colours in config: ",
        paste(
            missing_species_colors,
            collapse = ", "
        )
    )
}

if (
    length(padj_threshold) != 1L ||
    !is.finite(padj_threshold)
) {
    stop(
        "Invalid deseq2.padj_threshold in config."
    )
}

if (
    length(log2fc_threshold) != 1L ||
    !is.finite(log2fc_threshold)
) {
    stop(
        "Invalid deseq2.log2fc_threshold in config."
    )
}

if (
    length(show_category) != 1L ||
    !is.finite(show_category) ||
    show_category < 1
) {
    stop(
        "Invalid plots.general.showCategory in config."
    )
}

message(
    "DEG thresholds: adjusted p-value < ",
    padj_threshold,
    "; absolute log2 fold change >= ",
    log2fc_threshold
)

# ============================================================
# Validate inputs
# ============================================================

required_result_columns <- c(
    "GeneID",
    "padj",
    "log2FoldChange"
)

validate_results <- function(x, species_name) {

    missing_columns <- setdiff(
        required_result_columns,
        colnames(x)
    )

    if (length(missing_columns)) {
        stop(
            "Missing columns in ",
            species_name,
            ": ",
            paste(missing_columns, collapse = ", ")
        )
    }
}

validate_results(
    annot_exp,
    "experimentalis"
)

validate_results(
    annot_gad,
    "gadabouti"
)

required_term2gene <- c(
    "GO",
    "Gene"
)

required_term2name <- c(
    "GO",
    "Description"
)

if (
    any(
        !required_term2gene %in%
            colnames(term2gene)
    )
) {
    stop(
        "TERM2GENE must contain GO and Gene columns."
    )
}

if (
    any(
        !required_term2name %in%
            colnames(term2name)
    )
) {
    stop(
        "TERM2NAME must contain GO and Description columns."
    )
}

# ============================================================
# Clean annotation tables
# ============================================================

# LOC is preserved because both DESeq2 and TERM2GENE already use
# the same complete identifiers.

term2gene <- term2gene %>%
    dplyr::transmute(
        GO = trimws(
            as.character(GO)
        ),
        Gene = trimws(
            as.character(Gene)
        )
    ) %>%
    dplyr::filter(
        !is.na(GO),
        !is.na(Gene),
        GO != "",
        Gene != ""
    ) %>%
    dplyr::distinct()

term2name <- term2name %>%
    dplyr::transmute(
        GO = trimws(
            as.character(GO)
        ),
        Description = trimws(
            as.character(Description)
        )
    ) %>%
    dplyr::filter(
        !is.na(GO),
        !is.na(Description),
        GO != "",
        Description != ""
    ) %>%
    dplyr::distinct(
        GO,
        .keep_all = TRUE
    )

annot_exp$GeneID <- trimws(
    as.character(annot_exp$GeneID)
)

annot_gad$GeneID <- trimws(
    as.character(annot_gad$GeneID)
)

# ============================================================
# Classify genes using the final DEA thresholds
# ============================================================

classify_regulation <- function(x) {

    x %>%
        dplyr::mutate(
            regulation = dplyr::case_when(
                !is.na(padj) &
                    !is.na(log2FoldChange) &
                    padj < padj_threshold &
                    log2FoldChange >=
                        log2fc_threshold ~
                    "up",

                !is.na(padj) &
                    !is.na(log2FoldChange) &
                    padj < padj_threshold &
                    log2FoldChange <=
                        -log2fc_threshold ~
                    "down",

                TRUE ~ "nonsig"
            )
        )
}

annot_exp <- classify_regulation(
    annot_exp
)

annot_gad <- classify_regulation(
    annot_gad
)

# ============================================================
# Define gene sets and annotated backgrounds
# ============================================================

annotation_genes <- unique(
    term2gene$Gene
)

exp_background <- intersect(
    unique(annot_exp$GeneID),
    annotation_genes
)

gad_background <- intersect(
    unique(annot_gad$GeneID),
    annotation_genes
)

gene_sets <- list(
    experimentalis_up = intersect(
        unique(
            annot_exp$GeneID[
                annot_exp$regulation == "up"
            ]
        ),
        exp_background
    ),

    experimentalis_down = intersect(
        unique(
            annot_exp$GeneID[
                annot_exp$regulation == "down"
            ]
        ),
        exp_background
    ),

    gadabouti_up = intersect(
        unique(
            annot_gad$GeneID[
                annot_gad$regulation == "up"
            ]
        ),
        gad_background
    ),

    gadabouti_down = intersect(
        unique(
            annot_gad$GeneID[
                annot_gad$regulation == "down"
            ]
        ),
        gad_background
    )
)

backgrounds <- list(
    experimentalis_up = exp_background,
    experimentalis_down = exp_background,
    gadabouti_up = gad_background,
    gadabouti_down = gad_background
)

display_titles <- c(
    experimentalis_up =
        "Pam. experimentalis - Upregulated",

    experimentalis_down =
        "Pam. experimentalis - Downregulated",

    gadabouti_up =
        "Pam. gadabouti - Upregulated",

    gadabouti_down =
        "Pam. gadabouti - Downregulated"
)

message("")
message("Genes used for GO enrichment:")

for (analysis_name in names(gene_sets)) {
    message(
        analysis_name,
        ": ",
        length(gene_sets[[analysis_name]]),
        " genes; background = ",
        length(backgrounds[[analysis_name]])
    )
}

# ============================================================
# Run enrichment and save complete results
# ============================================================

enrichment_objects <- list()
enrichment_tables <- list()

for (analysis_name in names(gene_sets)) {

    genes <- gene_sets[[analysis_name]]
    universe <- backgrounds[[analysis_name]]

    message("")
    message(
        "Running enrichment: ",
        analysis_name
    )

    if (!length(genes)) {

        warning(
            "No genes available for ",
            analysis_name
        )

        enrichment_objects[[analysis_name]] <- NULL
        enrichment_tables[[analysis_name]] <-
            data.frame()

        saveRDS(
            NULL,
            rds_outputs[[analysis_name]]
        )

        utils::write.csv(
            data.frame(),
            csv_outputs[[analysis_name]],
            row.names = FALSE
        )

        next
    }

    ego <- tryCatch(
        clusterProfiler::enricher(
            gene = genes,
            universe = universe,
            TERM2GENE = term2gene,
            TERM2NAME = term2name,
            pAdjustMethod = "BH",
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.20,
            minGSSize = 10,
            maxGSSize = 500
        ),
        error = function(e) {
            warning(
                "Enricher failed for ",
                analysis_name,
                ": ",
                conditionMessage(e)
            )
            NULL
        }
    )

    if (
        is.null(ego) ||
        nrow(as.data.frame(ego)) == 0
    ) {

        warning(
            "No enrichment results for ",
            analysis_name
        )

        enrichment_objects[[analysis_name]] <- NULL
        enrichment_tables[[analysis_name]] <-
            data.frame()

        saveRDS(
            NULL,
            rds_outputs[[analysis_name]]
        )

        utils::write.csv(
            data.frame(),
            csv_outputs[[analysis_name]],
            row.names = FALSE
        )

        next
    }

    result_table <- as.data.frame(
        ego
    )

    result_table <- result_table %>%
        dplyr::mutate(
            minus_log10_padj = -log10(
                pmax(
                    p.adjust,
                    .Machine$double.xmin
                )
            )
        )

    enrichment_objects[[analysis_name]] <-
        ego

    enrichment_tables[[analysis_name]] <-
        result_table

    saveRDS(
        ego,
        rds_outputs[[analysis_name]]
    )

    utils::write.csv(
        result_table,
        csv_outputs[[analysis_name]],
        row.names = FALSE
    )

    message(
        analysis_name,
        ": ",
        nrow(result_table),
        " enriched GO terms saved"
    )
}

# ============================================================
# Determine a common -log10 adjusted p-value colour scale
# ============================================================

all_log_padj <- unlist(
    lapply(
        enrichment_tables,
        function(x) {

            if (
                is.null(x) ||
                !nrow(x) ||
                !"minus_log10_padj" %in%
                    colnames(x)
            ) {
                return(numeric())
            }

            x$minus_log10_padj[
                is.finite(
                    x$minus_log10_padj
                )
            ]
        }
    ),
    use.names = FALSE
)

if (!length(all_log_padj)) {
    colour_limits <- c(0, 1)
} else {
    colour_limits <- range(
        all_log_padj,
        na.rm = TRUE
    )

    if (
        !is.finite(colour_limits[1]) ||
        !is.finite(colour_limits[2]) ||
        colour_limits[1] ==
            colour_limits[2]
    ) {
        colour_limits <- c(
            0,
            max(
                1,
                colour_limits[2]
            )
        )
    }
}

# ============================================================
# Helpers for plotting
# ============================================================

parse_gene_ratio <- function(x) {

    vapply(
        strsplit(
            as.character(x),
            "/",
            fixed = TRUE
        ),
        function(parts) {

            if (length(parts) != 2L) {
                return(NA_real_)
            }

            numerator <- suppressWarnings(
                as.numeric(parts[1])
            )

            denominator <- suppressWarnings(
                as.numeric(parts[2])
            )

            if (
                !is.finite(numerator) ||
                !is.finite(denominator) ||
                denominator == 0
            ) {
                return(NA_real_)
            }

            numerator / denominator
        },
        numeric(1)
    )
}

make_empty_plot <- function(title_text) {

    ggplot() +
        annotate(
            "text",
            x = 0.5,
            y = 0.5,
            label = "No significant enrichment results",
            size = 5
        ) +
        xlim(0, 1) +
        ylim(0, 1) +
        labs(
            title = title_text
        ) +
        theme_void(
            base_size = 12
        ) +
        theme(
            plot.title = element_text(
                face = "bold",
                hjust = 0.5
            )
        )
}

make_bubble_plot <- function(
    result_table,
    title_text,
    species_color
) {

    if (
        is.null(result_table) ||
        !nrow(result_table)
    ) {
        return(
            make_empty_plot(
                title_text
            )
        )
    }

    plot_data <- result_table %>%
        dplyr::filter(
            !is.na(Description),
            !is.na(p.adjust),
            !is.na(Count),
            !is.na(GeneRatio)
        ) %>%
        dplyr::arrange(
            p.adjust,
            dplyr::desc(Count)
        ) %>%
        dplyr::slice_head(
            n = show_category
        ) %>%
        dplyr::arrange(
            dplyr::desc(Count),
            p.adjust
        ) %>%
        dplyr::mutate(
            gene_ratio_numeric =
                parse_gene_ratio(
                    GeneRatio
                ),

            minus_log10_padj = -log10(
                pmax(
                    p.adjust,
                    .Machine$double.xmin
                )
            ),

            Description = factor(
                Description,
                levels = rev(
                    unique(Description)
                )
            )
        ) %>%
        dplyr::filter(
            is.finite(
                gene_ratio_numeric
            ),
            is.finite(
                minus_log10_padj
            )
        )

    if (!nrow(plot_data)) {
        return(
            make_empty_plot(
                title_text
            )
        )
    }

    light_species_color <- scales::colour_ramp(
        c(
            "white",
            species_color
        )
    )(0.30)

    ggplot(
        plot_data,
        aes(
            x = gene_ratio_numeric,
            y = Description,
            size = Count,
            color = minus_log10_padj
        )
    ) +
        geom_point(
            alpha = 0.9
        ) +
        scale_color_gradient(
            low = light_species_color,
            high = species_color,
            limits = colour_limits,
            oob = scales::squish,
            name = expression(
                -log[10]("adjusted p-value")
            )
        ) +
        scale_size_continuous(
            name = "Gene count",
            range = c(3, 10)
        ) +
        scale_x_continuous(
            expand = expansion(
                mult = c(0.04, 0.22)
            )
        ) +
        coord_cartesian(
            clip = "off"
        ) +
        labs(
            title = title_text,
            x = "Gene ratio",
            y = NULL
        ) +
        guides(
            size = guide_legend(
                order = 1,
                override.aes = list(
                    color = species_color,
                    alpha = 0.9
                )
            ),
            color = guide_colorbar(
                order = 2
            )
        ) +
        theme_minimal(
            base_size = 12
        ) +
        theme(
            plot.title = element_text(
                face = "bold",
                hjust = 0.5,
                size = 14
            ),
            axis.text.y = element_text(
                size = 10
            ),
            panel.grid.major.y =
                element_blank(),
            panel.grid.minor =
                element_blank(),
            legend.position = "right",
            plot.margin = margin(
                8,
                30,
                8,
                8
            )
            
        )
}

# ============================================================
# Generate PDFs
# ============================================================

for (analysis_name in names(pdf_outputs)) {

    species_name <- sub(
        "_.*$",
        "",
        analysis_name
    )

    species_color <- species_colors[[species_name]]

    plot_object <- make_bubble_plot(
        result_table =
            enrichment_tables[[analysis_name]],

        title_text =
            display_titles[[analysis_name]],

        species_color =
            species_color
    )

    ggsave(
        filename = pdf_outputs[[analysis_name]],
        plot = plot_object,
        width = max(plot_width, 10),
        height = plot_height,
        device = cairo_pdf
    )

    message(
        "Saved plot: ",
        pdf_outputs[[analysis_name]]
    )
}

message("")
message("GO enrichment completed successfully.")