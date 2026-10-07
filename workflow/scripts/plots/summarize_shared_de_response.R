#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(readr)
})

# ============================================================
# Configuration
# ============================================================

required_species <- c(
    "experimentalis",
    "gadabouti"
)

category_order <- c(
    "Shared upregulated",
    "Shared downregulated",
    "Shared discordant",
    "Experimentalis-specific",
    "Gadabouti-specific"
)

response_colors <- c(
    "Shared upregulated" = "#009E73",
    "Shared downregulated" = "#0072B2",
    "Shared discordant" = "#D55E00",
    "Experimentalis-specific" = "#56b4e9",
    "Gadabouti-specific" = "#4b2588ff"
)

# ============================================================
# Read input
# ============================================================

res_list <- readRDS(
    snakemake@input[["deseq_results"]]
)

missing_species <- setdiff(
    required_species,
    names(res_list)
)

if (length(missing_species)) {
    stop(
        "Species missing from DESeq2 results: ",
        paste(missing_species, collapse = ", ")
    )
}

# ============================================================
# Prepare significant genes
# ============================================================

prepare_degs <- function(results, species_name) {

    required_columns <- c(
        "GeneID",
        "log2FoldChange",
        "sig"
    )

    missing_columns <- setdiff(
        required_columns,
        colnames(results)
    )

    if (length(missing_columns)) {
        stop(
            "Missing columns for ",
            species_name,
            ": ",
            paste(missing_columns, collapse = ", ")
        )
    }

    # Preserve annotation status when available
    if (!"annotation_status" %in% colnames(results)) {
        results$annotation_status <- NA_character_
    }

    results %>%
        filter(
            !is.na(GeneID),
            !is.na(log2FoldChange),
            !is.na(sig),
            sig == "yes"
        ) %>%
        distinct(
            GeneID,
            .keep_all = TRUE
        ) %>%
        transmute(
            GeneID = as.character(GeneID),
            log2FoldChange = as.numeric(log2FoldChange),
            annotation_status = as.character(annotation_status)
        )
}

experimentalis <- prepare_degs(
    res_list[["experimentalis"]],
    "experimentalis"
) %>%
    rename(
        log2FC_experimentalis = log2FoldChange,
        annotation_experimentalis = annotation_status
    )

gadabouti <- prepare_degs(
    res_list[["gadabouti"]],
    "gadabouti"
) %>%
    rename(
        log2FC_gadabouti = log2FoldChange,
        annotation_gadabouti = annotation_status
    )

n_experimentalis <- nrow(experimentalis)
n_gadabouti <- nrow(gadabouti)

message(
    "Experimentalis DEGs: ",
    n_experimentalis
)

message(
    "Gadabouti DEGs: ",
    n_gadabouti
)

# ============================================================
# Classify shared and species-specific responses
# ============================================================

gene_classification <- full_join(
    experimentalis,
    gadabouti,
    by = "GeneID"
) %>%
    mutate(
        significant_experimentalis = !is.na(
            log2FC_experimentalis
        ),
        significant_gadabouti = !is.na(
            log2FC_gadabouti
        ),

        response_category = case_when(
            significant_experimentalis &
                !significant_gadabouti ~
                "Experimentalis-specific",

            !significant_experimentalis &
                significant_gadabouti ~
                "Gadabouti-specific",

            significant_experimentalis &
                significant_gadabouti &
                log2FC_experimentalis > 0 &
                log2FC_gadabouti > 0 ~
                "Shared upregulated",

            significant_experimentalis &
                significant_gadabouti &
                log2FC_experimentalis < 0 &
                log2FC_gadabouti < 0 ~
                "Shared downregulated",

            significant_experimentalis &
                significant_gadabouti &
                log2FC_experimentalis *
                log2FC_gadabouti < 0 ~
                "Shared discordant",

            TRUE ~ "Unclassified"
        ),

        shared_response = case_when(
            response_category == "Shared upregulated" ~
                "Concordant",

            response_category == "Shared downregulated" ~
                "Concordant",

            response_category == "Shared discordant" ~
                "Discordant",

            TRUE ~ "Species-specific"
        ),

        annotation_status = coalesce(
            annotation_experimentalis,
            annotation_gadabouti
        )
    ) %>%
    arrange(
        factor(
            response_category,
            levels = c(
                category_order,
                "Unclassified"
            )
        ),
        GeneID
    )

unclassified <- gene_classification %>%
    filter(
        response_category == "Unclassified"
    )

if (nrow(unclassified) > 0) {
    warning(
        nrow(unclassified),
        " significant genes could not be classified. ",
        "This usually occurs when a significant gene has a ",
        "log2 fold change equal to zero."
    )
}

# ============================================================
# Response-category summary
# ============================================================

classified_genes <- gene_classification %>%
    filter(
        response_category != "Unclassified"
    )

n_union <- nrow(classified_genes)

response_summary <- classified_genes %>%
    count(
        response_category,
        name = "number_of_genes"
    ) %>%
    mutate(
        response_category = factor(
            response_category,
            levels = category_order
        ),
        percentage_of_de_union = round(
            100 * number_of_genes / n_union,
            2
        )
    ) %>%
    arrange(response_category) %>%
    mutate(
        response_category = as.character(
            response_category
        )
    )

# Add absent categories with zero genes
response_summary <- tibble(
    response_category = category_order
) %>%
    left_join(
        response_summary,
        by = "response_category"
    ) %>%
    mutate(
        number_of_genes = coalesce(
            number_of_genes,
            0L
        ),
        percentage_of_de_union = coalesce(
            percentage_of_de_union,
            0
        )
    )

# ============================================================
# General overlap summary
# ============================================================

n_shared <- sum(
    gene_classification$significant_experimentalis &
        gene_classification$significant_gadabouti
)

n_experimentalis_specific <- sum(
    gene_classification$response_category ==
        "Experimentalis-specific"
)

n_gadabouti_specific <- sum(
    gene_classification$response_category ==
        "Gadabouti-specific"
)

n_concordant_up <- sum(
    gene_classification$response_category ==
        "Shared upregulated"
)

n_concordant_down <- sum(
    gene_classification$response_category ==
        "Shared downregulated"
)

n_discordant <- sum(
    gene_classification$response_category ==
        "Shared discordant"
)

overlap_summary <- tibble(
    category = c(
        "Experimentalis DEGs",
        "Gadabouti DEGs",
        "Shared DEGs",
        "Experimentalis-specific DEGs",
        "Gadabouti-specific DEGs",
        "Shared concordant upregulated",
        "Shared concordant downregulated",
        "Shared discordant"
    ),
    number_of_genes = c(
        n_experimentalis,
        n_gadabouti,
        n_shared,
        n_experimentalis_specific,
        n_gadabouti_specific,
        n_concordant_up,
        n_concordant_down,
        n_discordant
    )
)

message(
    "Shared DEGs: ",
    n_shared
)

message(
    "Experimentalis-specific DEGs: ",
    n_experimentalis_specific
)

message(
    "Gadabouti-specific DEGs: ",
    n_gadabouti_specific
)

message(
    "Shared concordant upregulated: ",
    n_concordant_up
)

message(
    "Shared concordant downregulated: ",
    n_concordant_down
)

message(
    "Shared discordant: ",
    n_discordant
)

# ============================================================
# Create output directories
# ============================================================

output_paths <- c(
    as.character(
        snakemake@output[["gene_table"]]
    ),
    as.character(
        snakemake@output[["response_summary"]]
    ),
    as.character(
        snakemake@output[["overlap_summary"]]
    ),
    as.character(
        snakemake@output[["plot"]]
    )
)

invisible(
    lapply(
        unique(dirname(output_paths)),
        dir.create,
        recursive = TRUE,
        showWarnings = FALSE
    )
)

# ============================================================
# Write output tables
# ============================================================

write_csv(
    gene_classification,
    snakemake@output[["gene_table"]],
    na = ""
)

write_csv(
    response_summary,
    snakemake@output[["response_summary"]],
    na = ""
)

write_csv(
    overlap_summary,
    snakemake@output[["overlap_summary"]],
    na = ""
)

# ============================================================
# Plot response categories
# ============================================================

plot_data <- response_summary %>%
    mutate(
        response_category = factor(
            response_category,
            levels = rev(category_order)
        ),
        label = paste0(
            number_of_genes,
            " (",
            percentage_of_de_union,
            "%)"
        )
    )

p <- ggplot(
    plot_data,
    aes(
        x = number_of_genes,
        y = response_category,
        fill = response_category
    )
) +
    geom_col(
        width = 0.72
    ) +
    geom_text(
        aes(label = label),
        hjust = -0.1,
        size = 3.5
    ) +
    scale_fill_manual(
        values = response_colors,
        guide = "none"
    ) +
    scale_x_continuous(
        expand = expansion(
            mult = c(0, 0.25)
        )
    ) +
    labs(
        x = "Number of differentially expressed genes",
        y = NULL
    ) +
    theme_minimal(
        base_size = 11
    ) +
    theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(
            size = 10
        ),
        axis.title.x = element_text(
            size = 10
        ),
        plot.margin = margin(
            10,
            20,
            10,
            10
        )
    )

ggsave(
    filename = snakemake@output[["plot"]],
    plot = p,
    width = 8,
    height = 5,
    device = cairo_pdf
)