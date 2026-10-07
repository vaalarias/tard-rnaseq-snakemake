#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(DESeq2)
    library(dplyr)
    library(tibble)
    library(yaml)
})

# ============================================================
# 1. Read configuration
# ============================================================

config_path <- as.character(
    snakemake@input[["config"]]
)

raw_counts_path <- as.character(
    snakemake@input[["raw_counts"]]
)

config <- yaml::read_yaml(
    config_path
)

design_formula <- as.formula(
    config$deseq2$design
)

contrast_vector <- as.character(
    config$deseq2$contrasts[[1]]
)

padj_threshold <- as.numeric(
    config$deseq2$padj_threshold
)

log2fc_threshold <- as.numeric(
    config$deseq2$log2fc_threshold
)

if (length(contrast_vector) != 3) {
    stop(
        "The DESeq2 contrast must contain exactly three elements: ",
        "factor, numerator level, and denominator level. Found: ",
        paste(contrast_vector, collapse = ", ")
    )
}

contrast_factor <- contrast_vector[1]
contrast_numerator <- contrast_vector[2]
contrast_denominator <- contrast_vector[3]

message(
    "Config file: ",
    config_path
)

message(
    "DESeq2 design: ",
    deparse(design_formula)
)

message(
    "DESeq2 contrast: ",
    contrast_factor,
    ", ",
    contrast_numerator,
    " versus ",
    contrast_denominator
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
# 2. Read count matrix and sample metadata
# ============================================================

raw <- readRDS(
    raw_counts_path
)

cts <- as.matrix(
    raw$counts
)

coldata <- as.data.frame(
    raw$coldata
)

required_metadata <- c(
    "species",
    contrast_factor
)

missing_metadata <- setdiff(
    required_metadata,
    colnames(coldata)
)

if (length(missing_metadata)) {
    stop(
        "Missing columns in coldata: ",
        paste(missing_metadata, collapse = ", ")
    )
}

if (is.null(rownames(coldata))) {
    stop(
        "coldata must have sample identifiers as row names."
    )
}

missing_count_samples <- setdiff(
    colnames(cts),
    rownames(coldata)
)

missing_metadata_samples <- setdiff(
    rownames(coldata),
    colnames(cts)
)

if (length(missing_count_samples)) {
    stop(
        "Count-matrix samples missing from coldata: ",
        paste(missing_count_samples, collapse = ", ")
    )
}

if (length(missing_metadata_samples)) {
    warning(
        "Metadata samples absent from the count matrix: ",
        paste(missing_metadata_samples, collapse = ", ")
    )
}

# Put coldata in exactly the same order as the count matrix
coldata <- coldata[
    colnames(cts),
    ,
    drop = FALSE
]

# Explicitly convert variables to factors
coldata$species <- factor(
    coldata$species
)

coldata[[contrast_factor]] <- factor(
    coldata[[contrast_factor]]
)

available_conditions <- levels(
    coldata[[contrast_factor]]
)

missing_contrast_levels <- setdiff(
    c(
        contrast_numerator,
        contrast_denominator
    ),
    available_conditions
)

if (length(missing_contrast_levels)) {
    stop(
        "Contrast levels missing from coldata: ",
        paste(missing_contrast_levels, collapse = ", "),
        ". Available levels: ",
        paste(available_conditions, collapse = ", ")
    )
}

# Set control as reference for clearer coefficient names.
# The explicit contrast below still determines the reported direction.
coldata[[contrast_factor]] <- relevel(
    coldata[[contrast_factor]],
    ref = contrast_denominator
)

# ============================================================
# 3. Run DESeq2 independently for each species
# ============================================================

species_list <- levels(
    droplevels(coldata$species)
)

results_all <- setNames(
    vector(
        mode = "list",
        length = length(species_list)
    ),
    species_list
)

dds_all <- setNames(
    vector(
        mode = "list",
        length = length(species_list)
    ),
    species_list
)

for (sp in species_list) {

    message(
        "\n========================================"
    )

    message(
        "Running DESeq2 for species: ",
        sp
    )

    sample_ids <- rownames(coldata)[
        coldata$species == sp
    ]

    if (!length(sample_ids)) {
        stop(
            "No samples found for species: ",
            sp
        )
    }

    counts_sp <- cts[
        ,
        sample_ids,
        drop = FALSE
    ]

    coldata_sp <- droplevels(
        coldata[
            sample_ids,
            ,
            drop = FALSE
        ]
    )

    message(
        "Samples: ",
        paste(sample_ids, collapse = ", ")
    )

    message(
        "Condition table:"
    )

    print(
        table(
            coldata_sp[[contrast_factor]]
        )
    )

    species_levels <- levels(
        coldata_sp[[contrast_factor]]
    )

    missing_species_levels <- setdiff(
        c(
            contrast_numerator,
            contrast_denominator
        ),
        species_levels
    )

    if (length(missing_species_levels)) {
        stop(
            "Species ",
            sp,
            " does not contain all contrast levels. Missing: ",
            paste(missing_species_levels, collapse = ", ")
        )
    }

    # Construct the species-specific DESeq2 object directly
    dds_sp <- DESeqDataSetFromMatrix(
        countData = counts_sp,
        colData = coldata_sp,
        design = design_formula
    )

    # No count prefilter is applied here. DESeq2 handles low counts
    # through independent filtering when results() is called.
    dds_sp <- DESeq(
        dds_sp
    )

    message(
        "Available DESeq2 coefficients: ",
        paste(
            resultsNames(dds_sp),
            collapse = ", "
        )
    )

    res <- results(
        dds_sp,
        contrast = contrast_vector
    )

    res_df <- as.data.frame(
        res
    ) %>%
        rownames_to_column(
            "GeneID"
        ) %>%
        mutate(
            sig = case_when(
                !is.na(padj) &
                    padj < padj_threshold &
                    abs(log2FoldChange) >=
                    log2fc_threshold ~
                    "yes",

                TRUE ~ "no"
            ),
            direction = case_when(
                sig == "yes" &
                    log2FoldChange > 0 ~
                    "upregulated",

                sig == "yes" &
                    log2FoldChange < 0 ~
                    "downregulated",

                TRUE ~ "not_significant"
            )
        )

    n_tested <- sum(
        !is.na(res_df$padj)
    )

    n_up <- sum(
        res_df$direction == "upregulated"
    )

    n_down <- sum(
        res_df$direction == "downregulated"
    )

    message(
        "Genes with non-NA adjusted p-values: ",
        n_tested
    )

    message(
        "Upregulated DEGs: ",
        n_up
    )

    message(
        "Downregulated DEGs: ",
        n_down
    )

    message(
        "Total DEGs: ",
        n_up + n_down
    )

    results_all[[sp]] <- res_df
    dds_all[[sp]] <- dds_sp
}

# ============================================================
# 4. Save results
# ============================================================

deseq_output <- as.character(
    snakemake@output[["deseq_results"]]
)

dir.create(
    dirname(deseq_output),
    recursive = TRUE,
    showWarnings = FALSE
)

saveRDS(
    results_all,
    deseq_output
)

# Optional: save fitted species-specific DESeq2 objects if the
# rule defines a named output called "dds_list".
if ("dds_list" %in% names(snakemake@output)) {

    dds_output <- as.character(
        snakemake@output[["dds_list"]]
    )

    dir.create(
        dirname(dds_output),
        recursive = TRUE,
        showWarnings = FALSE
    )

    saveRDS(
        dds_all,
        dds_output
    )
}

# ============================================================
# 5. Save species-specific CSV files
# ============================================================

out_dir <- as.character(
    snakemake@params[["outdir"]]
)

dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

for (sp in names(results_all)) {

    csv_path <- file.path(
        out_dir,
        paste0(
            sp,
            "_deseq2_results.csv"
        )
    )

    write.csv(
        results_all[[sp]],
        csv_path,
        row.names = FALSE
    )

    message(
        "Written: ",
        csv_path
    )
}