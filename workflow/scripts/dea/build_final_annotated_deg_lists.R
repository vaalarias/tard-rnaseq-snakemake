#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(readr)
    library(yaml)
})

# ============================================================
# Inputs
# ============================================================

exp_results <- read.csv(
    snakemake@input[["exp_results"]],
    stringsAsFactors = FALSE,
    check.names = FALSE
)

gad_results <- read.csv(
    snakemake@input[["gad_results"]],
    stringsAsFactors = FALSE,
    check.names = FALSE
)

gaf_file <- snakemake@input[["gaf"]]
gtf_file <- snakemake@input[["gtf"]]

term2name <- readRDS(
    snakemake@input[["term2name"]]
)

config <- yaml::read_yaml(
    snakemake@input[["config"]]
)

exp_output <- snakemake@output[["exp"]]
gad_output <- snakemake@output[["gad"]]

dir.create(
    dirname(exp_output),
    recursive = TRUE,
    showWarnings = FALSE
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

# ============================================================
# Helpers
# ============================================================

clean_gene_id <- function(x) {
    trimws(
        as.character(x)
    )
}

collapse_unique <- function(x) {

    x <- trimws(
        as.character(x)
    )

    x <- x[
        !is.na(x) &
        nzchar(x) &
        x != "-"
    ]

    paste(
        unique(x),
        collapse = ";"
    )
}

extract_gtf_attribute <- function(
    attributes,
    attribute_name
) {

    pattern <- paste0(
        "(?:^|;[[:space:]]*)",
        attribute_name,
        "[[:space:]]+\"([^\"]*)\""
    )

    match <- stringr::str_match(
        attributes,
        pattern
    )

    match[, 2]
}

select_degs <- function(results, species_name) {

    required_columns <- c(
        "GeneID",
        "log2FoldChange",
        "padj"
    )

    missing_columns <- setdiff(
        required_columns,
        colnames(results)
    )

    if (length(missing_columns)) {
        stop(
            "Missing columns in ",
            species_name,
            ": ",
            paste(
                missing_columns,
                collapse = ", "
            )
        )
    }

    if ("sig" %in% colnames(results)) {

        degs <- results %>%
            filter(
                !is.na(sig),
                sig == "yes"
            )

    } else {

        degs <- results %>%
            filter(
                !is.na(padj),
                !is.na(log2FoldChange),
                padj < padj_threshold,
                abs(log2FoldChange) >=
                    log2fc_threshold
            )
    }

    degs %>%
        transmute(
            Gene_Symbol = clean_gene_id(
                GeneID
            ),
            log2FoldChange,
            padj
        ) %>%
        distinct(
            Gene_Symbol,
            .keep_all = TRUE
        ) %>%
        arrange(
            desc(log2FoldChange)
        )
}

# ============================================================
# Read direct GAF annotations
# ============================================================

gaf_columns <- c(
    "DB",
    "DB_Object_ID",
    "DB_Object_Symbol",
    "Qualifier",
    "GO_ID",
    "DB_Reference",
    "Evidence_Code",
    "With_From",
    "Aspect",
    "DB_Object_Name",
    "DB_Object_Synonym",
    "DB_Object_Type",
    "Taxon",
    "Date",
    "Assigned_By",
    "Annotation_Extension",
    "Gene_Product_Form_ID"
)

gaf <- read.delim(
    gaf_file,
    comment.char = "!",
    header = FALSE,
    sep = "\t",
    quote = "",
    fill = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (ncol(gaf) > length(gaf_columns)) {
    stop(
        "The GAF contains more than 17 columns."
    )
}

colnames(gaf) <- gaf_columns[
    seq_len(ncol(gaf))
]

# Ignore explicitly negated GO annotations.
gaf <- gaf %>%
    filter(
        is.na(Qualifier) |
        !str_detect(
            Qualifier,
            "(^|\\|)NOT($|\\|)"
        )
    )

# Allow either GAF object ID or symbol to match the LOC identifier.
gaf_long <- bind_rows(
    gaf %>%
        transmute(
            Gene_Symbol = clean_gene_id(
                DB_Object_ID
            ),
            GO_ID = as.character(GO_ID)
        ),

    gaf %>%
        transmute(
            Gene_Symbol = clean_gene_id(
                DB_Object_Symbol
            ),
            GO_ID = as.character(GO_ID)
        )
) %>%
    filter(
        !is.na(Gene_Symbol),
        nzchar(Gene_Symbol),
        !is.na(GO_ID),
        nzchar(GO_ID)
    ) %>%
    distinct()

# ============================================================
# GO names
# ============================================================

term2name <- as.data.frame(
    term2name,
    stringsAsFactors = FALSE
)

if (
    all(
        c("GO", "Description") %in%
            colnames(term2name)
    )
) {

    go_names <- term2name %>%
        transmute(
            GO_ID = as.character(GO),
            GO_Term_Description =
                as.character(Description)
        )

} else {

    go_names <- term2name[, 1:2] %>%
        setNames(
            c(
                "GO_ID",
                "GO_Term_Description"
            )
        )
}

go_names <- go_names %>%
    filter(
        !is.na(GO_ID),
        nzchar(GO_ID)
    ) %>%
    distinct(
        GO_ID,
        .keep_all = TRUE
    )

gene_go_annotations <- gaf_long %>%
    left_join(
        go_names,
        by = "GO_ID"
    ) %>%
    group_by(
        Gene_Symbol
    ) %>%
    summarise(
        GO_Terms = collapse_unique(
            GO_ID
        ),
        GO_Term_Description =
            collapse_unique(
                GO_Term_Description
            ),
        .groups = "drop"
    )

# ============================================================
# Read GTF annotations
# ============================================================

gtf <- read.delim(
    gtf_file,
    comment.char = "#",
    header = FALSE,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE,
    fill = TRUE
)

if (ncol(gtf) < 9) {
    stop(
        "Invalid GTF: fewer than nine columns."
    )
}

colnames(gtf)[1:9] <- c(
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attributes"
)

gtf_annotations <- gtf %>%
    transmute(
        Gene_Symbol = clean_gene_id(
            extract_gtf_attribute(
                attributes,
                "gene_id"
            )
        ),

        gene_name = extract_gtf_attribute(
            attributes,
            "gene"
        ),

        product = extract_gtf_attribute(
            attributes,
            "product"
        ),

        description =
            extract_gtf_attribute(
                attributes,
                "description"
            ),

        Protein_ID =
            extract_gtf_attribute(
                attributes,
                "protein_id"
            )
    ) %>%
    filter(
        !is.na(Gene_Symbol),
        nzchar(Gene_Symbol)
    ) %>%
    group_by(
        Gene_Symbol
    ) %>%
    summarise(
        GTF_annotation = collapse_unique(
            c(
                gene_name,
                product,
                description
            )
        ),

        Protein_ID = collapse_unique(
            Protein_ID
        ),

        .groups = "drop"
    )

# ============================================================
# Construct final tables
# ============================================================

build_final_table <- function(
    results,
    species_name
) {

    degs <- select_degs(
        results,
        species_name
    )

    final <- degs %>%
        left_join(
            gene_go_annotations,
            by = "Gene_Symbol"
        ) %>%
        left_join(
            gtf_annotations,
            by = "Gene_Symbol"
        ) %>%
        mutate(
            across(
                c(
                    GO_Terms,
                    GO_Term_Description,
                    GTF_annotation,
                    Protein_ID
                ),
                ~ replace_na(
                    as.character(.x),
                    ""
                )
            )
        ) %>%
        select(
            Gene_Symbol,
            log2FoldChange,
            padj,
            GO_Terms,
            GO_Term_Description,
            GTF_annotation,
            Protein_ID
        )

    if (
        nrow(final) !=
        n_distinct(final$Gene_Symbol)
    ) {
        stop(
            "Duplicated genes in final ",
            species_name,
            " table."
        )
    }

    message(
        species_name,
        ": ",
        nrow(final),
        " DEGs; ",
        sum(final$GO_Terms != ""),
        " with direct GO annotations; ",
        sum(final$GO_Terms == ""),
        " without direct GO annotations"
    )

    final
}

exp_final <- build_final_table(
    exp_results,
    "experimentalis"
)

gad_final <- build_final_table(
    gad_results,
    "gadabouti"
)

# ============================================================
# Validation
# ============================================================

if (nrow(exp_final) != 2146) {
    warning(
        "Expected 2146 experimentalis DEGs, obtained ",
        nrow(exp_final)
    )
}

if (nrow(gad_final) != 3911) {
    warning(
        "Expected 3911 gadabouti DEGs, obtained ",
        nrow(gad_final)
    )
}

# write.csv2 reproduces the previous semicolon separator and
# decimal-comma representation.
write.csv2(
    exp_final,
    exp_output,
    row.names = FALSE,
    quote = TRUE
)

write.csv2(
    gad_final,
    gad_output,
    row.names = FALSE,
    quote = TRUE
)