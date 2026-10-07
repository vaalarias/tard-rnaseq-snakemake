#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(yaml)
})

# ============================================================
# Inputs
# ============================================================

exp_results <- read.csv(
    snakemake@input[["exp_results"]],
    stringsAsFactors = FALSE
)

gad_results <- read.csv(
    snakemake@input[["gad_results"]],
    stringsAsFactors = FALSE
)

config <- yaml::read_yaml(
    snakemake@input[["config"]]
)

padj_threshold <- as.numeric(
    config$deseq2$padj_threshold
)

log2fc_threshold <- as.numeric(
    config$deseq2$log2fc_threshold
)

message(
    "Adjusted p-value threshold: ",
    padj_threshold
)

message(
    "Absolute log2 fold-change threshold: ",
    log2fc_threshold
)

# ============================================================
# Validate columns
# ============================================================

required_columns <- c(
    "GeneID",
    "padj",
    "log2FoldChange"
)

validate_results <- function(x, species) {

    missing_columns <- setdiff(
        required_columns,
        colnames(x)
    )

    if (length(missing_columns) > 0) {
        stop(
            "Missing columns in ",
            species,
            ": ",
            paste(
                missing_columns,
                collapse = ", "
            )
        )
    }

    if (anyDuplicated(x$GeneID)) {
        stop(
            "Duplicated GeneID values found in ",
            species,
            " results."
        )
    }
}

validate_results(
    exp_results,
    "experimentalis"
)

validate_results(
    gad_results,
    "gadabouti"
)

# ============================================================
# Classify differential expression
# ============================================================

classify_degs <- function(x) {

    x %>%
        mutate(
            GeneID = as.character(GeneID),

            regulation = case_when(
                !is.na(padj) &
                    padj < padj_threshold &
                    log2FoldChange >=
                    log2fc_threshold ~
                    "up",

                !is.na(padj) &
                    padj < padj_threshold &
                    log2FoldChange <=
                    -log2fc_threshold ~
                    "down",

                TRUE ~ "not_significant"
            )
        )
}

exp_results <- classify_degs(
    exp_results
)

gad_results <- classify_degs(
    gad_results
)

# ============================================================
# Generate gene sets
# ============================================================

gene_sets <- list(
    experimentalis = list(
        up = unique(
            exp_results$GeneID[
                !is.na(exp_results$GeneID) &
                    exp_results$regulation == "up"
            ]
        ),

        down = unique(
            exp_results$GeneID[
                !is.na(exp_results$GeneID) &
                    exp_results$regulation == "down"
            ]
        )
    ),

    gadabouti = list(
        up = unique(
            gad_results$GeneID[
                !is.na(gad_results$GeneID) &
                    gad_results$regulation == "up"
            ]
        ),

        down = unique(
            gad_results$GeneID[
                !is.na(gad_results$GeneID) &
                    gad_results$regulation == "down"
            ]
        )
    )
)

gene_sets$experimentalis$all <- union(
    gene_sets$experimentalis$up,
    gene_sets$experimentalis$down
)

gene_sets$gadabouti$all <- union(
    gene_sets$gadabouti$up,
    gene_sets$gadabouti$down
)

# ============================================================
# Shared genes
# ============================================================

gene_sets$shared <- list(
    up = intersect(
        gene_sets$experimentalis$up,
        gene_sets$gadabouti$up
    ),

    down = intersect(
        gene_sets$experimentalis$down,
        gene_sets$gadabouti$down
    ),

    all = intersect(
        gene_sets$experimentalis$all,
        gene_sets$gadabouti$all
    ),

    exp_up_gad_down = intersect(
        gene_sets$experimentalis$up,
        gene_sets$gadabouti$down
    ),

    exp_down_gad_up = intersect(
        gene_sets$experimentalis$down,
        gene_sets$gadabouti$up
    )
)

gene_sets$shared$same_direction <- union(
    gene_sets$shared$up,
    gene_sets$shared$down
)

gene_sets$shared$opposite_direction <- union(
    gene_sets$shared$exp_up_gad_down,
    gene_sets$shared$exp_down_gad_up
)

# ============================================================
# Calculate global totals
# ============================================================

n_experimentalis <- length(
    gene_sets$experimentalis$all
)

n_gadabouti <- length(
    gene_sets$gadabouti$all
)

n_shared <- length(
    gene_sets$shared$all
)

experimentalis_specific <- setdiff(
    gene_sets$experimentalis$all,
    gene_sets$gadabouti$all
)

gadabouti_specific <- setdiff(
    gene_sets$gadabouti$all,
    gene_sets$experimentalis$all
)

n_experimentalis_specific <- length(
    experimentalis_specific
)

n_gadabouti_specific <- length(
    gadabouti_specific
)

all_deg_union <- union(
    gene_sets$experimentalis$all,
    gene_sets$gadabouti$all
)

n_union <- length(
    all_deg_union
)

# ============================================================
# Validate overlap calculations
# ============================================================

if (
    n_experimentalis !=
        n_shared +
        n_experimentalis_specific
) {
    stop(
        "Experimentalis total is inconsistent with ",
        "shared + species-specific DEGs."
    )
}

if (
    n_gadabouti !=
        n_shared +
        n_gadabouti_specific
) {
    stop(
        "Gadabouti total is inconsistent with ",
        "shared + species-specific DEGs."
    )
}

if (
    n_union !=
        n_shared +
        n_experimentalis_specific +
        n_gadabouti_specific
) {
    stop(
        "The DEG union is inconsistent with the ",
        "overlap calculations."
    )
}

if (
    n_shared !=
        length(gene_sets$shared$same_direction) +
        length(gene_sets$shared$opposite_direction)
) {
    stop(
        "Shared DEGs are inconsistent with concordant + ",
        "discordant DEGs."
    )
}

# ============================================================
# Summary by species
# ============================================================

species_summary <- data.frame(
    species = c(
        "experimentalis",
        "gadabouti"
    ),

    upregulated = c(
        length(
            gene_sets$experimentalis$up
        ),
        length(
            gene_sets$gadabouti$up
        )
    ),

    downregulated = c(
        length(
            gene_sets$experimentalis$down
        ),
        length(
            gene_sets$gadabouti$down
        )
    ),

    total_degs = c(
        n_experimentalis,
        n_gadabouti
    ),

    shared_degs = c(
        n_shared,
        n_shared
    ),

    species_specific_degs = c(
        n_experimentalis_specific,
        n_gadabouti_specific
    )
) %>%
    mutate(
        percentage_shared = round(
            100 *
                shared_degs /
                total_degs,
            2
        ),

        percentage_species_specific = round(
            100 *
                species_specific_degs /
                total_degs,
            2
        )
    )

# ============================================================
# Shared and species-specific summary
# ============================================================

shared_summary <- data.frame(
    category = c(
        "Total experimentalis DEGs",
        "Total gadabouti DEGs",
        "DEGs in union",
        "Shared DEGs regardless of direction",
        "Experimentalis-specific DEGs",
        "Gadabouti-specific DEGs",
        "Shared upregulated",
        "Shared downregulated",
        "Shared same direction",
        "Shared opposite direction",
        "Experimentalis up / Gadabouti down",
        "Experimentalis down / Gadabouti up"
    ),

    number_of_genes = c(
        n_experimentalis,
        n_gadabouti,
        n_union,
        n_shared,
        n_experimentalis_specific,
        n_gadabouti_specific,

        length(
            gene_sets$shared$up
        ),

        length(
            gene_sets$shared$down
        ),

        length(
            gene_sets$shared$same_direction
        ),

        length(
            gene_sets$shared$opposite_direction
        ),

        length(
            gene_sets$shared$exp_up_gad_down
        ),

        length(
            gene_sets$shared$exp_down_gad_up
        )
    )
)

# ============================================================
# Gene-level comparison table
# ============================================================

all_deg_ids <- union(
    gene_sets$experimentalis$all,
    gene_sets$gadabouti$all
)

exp_comparison <- exp_results %>%
    select(
        GeneID,
        experimentalis_log2FC =
            log2FoldChange,
        experimentalis_padj =
            padj,
        experimentalis_regulation =
            regulation
    )

gad_comparison <- gad_results %>%
    select(
        GeneID,
        gadabouti_log2FC =
            log2FoldChange,
        gadabouti_padj =
            padj,
        gadabouti_regulation =
            regulation
    )

gene_comparison <- data.frame(
    GeneID = all_deg_ids,
    stringsAsFactors = FALSE
) %>%
    left_join(
        exp_comparison,
        by = "GeneID"
    ) %>%
    left_join(
        gad_comparison,
        by = "GeneID"
    ) %>%
    mutate(
        comparison_category = case_when(
            experimentalis_regulation == "up" &
                gadabouti_regulation == "up" ~
                "shared_upregulated",

            experimentalis_regulation == "down" &
                gadabouti_regulation == "down" ~
                "shared_downregulated",

            experimentalis_regulation == "up" &
                gadabouti_regulation == "down" ~
                "shared_discordant_exp_up_gad_down",

            experimentalis_regulation == "down" &
                gadabouti_regulation == "up" ~
                "shared_discordant_exp_down_gad_up",

            experimentalis_regulation %in%
                c("up", "down") &
                (
                    is.na(gadabouti_regulation) |
                    gadabouti_regulation ==
                        "not_significant"
                ) ~
                "experimentalis_specific",

            gadabouti_regulation %in%
                c("up", "down") &
                (
                    is.na(
                        experimentalis_regulation
                    ) |
                    experimentalis_regulation ==
                        "not_significant"
                ) ~
                "gadabouti_specific",

            TRUE ~ "unclassified"
        )
    ) %>%
    arrange(
        comparison_category,
        GeneID
    )

# ============================================================
# Validate gene-level table
# ============================================================

if (nrow(gene_comparison) != n_union) {
    stop(
        "The gene-level comparison table contains ",
        nrow(gene_comparison),
        " rows, but the expected DEG union contains ",
        n_union,
        " genes."
    )
}

n_unclassified <- sum(
    gene_comparison$comparison_category ==
        "unclassified"
)

if (n_unclassified > 0) {
    warning(
        n_unclassified,
        " genes could not be assigned to a comparison category."
    )
}

# ============================================================
# Create output directories
# ============================================================

output_paths <- c(
    as.character(
        snakemake@output[["species_summary"]]
    ),
    as.character(
        snakemake@output[["shared_summary"]]
    ),
    as.character(
        snakemake@output[["gene_comparison"]]
    ),
    as.character(
        snakemake@output[["gene_sets"]]
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
# Save outputs
# ============================================================

write.csv(
    species_summary,
    snakemake@output[["species_summary"]],
    row.names = FALSE
)

write.csv(
    shared_summary,
    snakemake@output[["shared_summary"]],
    row.names = FALSE
)

write.csv(
    gene_comparison,
    snakemake@output[["gene_comparison"]],
    row.names = FALSE
)

saveRDS(
    gene_sets,
    snakemake@output[["gene_sets"]]
)

# ============================================================
# Print summary
# ============================================================

message(
    "\nDifferential expression summary:"
)

print(
    species_summary
)

message(
    "\nShared DEG summary:"
)

print(
    shared_summary
)

message(
    "\nGene-level category counts:"
)

print(
    table(
        gene_comparison$comparison_category
    )
)
