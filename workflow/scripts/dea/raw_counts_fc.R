#!/usr/bin/env Rscript

library(tidyverse)

### ---------------------------
### 1. Read samplesheet
### ---------------------------

samples <- read.csv(
    snakemake@input[["samplesheet"]],
    comment.char = "#",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (!"sample" %in% colnames(samples)) {
    stop(
        "The samplesheet must contain a 'sample' column. Found: ",
        paste(colnames(samples), collapse = ", ")
    )
}

samples_unique <- samples[
    !duplicated(samples$sample),
    ,
    drop = FALSE
]

rownames(samples_unique) <- samples_unique$sample


### ---------------------------
### 2. Get featureCounts files
### ---------------------------

count_files <- unlist(
    snakemake@input[["counts"]],
    use.names = FALSE
)

if (length(count_files) == 0) {
    stop("No featureCounts files were provided by Snakemake.")
}

missing_files <- count_files[!file.exists(count_files)]

if (length(missing_files) > 0) {
    stop(
        "Missing featureCounts files: ",
        paste(missing_files, collapse = ", ")
    )
}

sample_names <- sub(
    "_counts\\.txt$",
    "",
    basename(count_files)
)

if (anyDuplicated(sample_names)) {
    stop(
        "Duplicated sample names obtained from count filenames: ",
        paste(sample_names[duplicated(sample_names)], collapse = ", ")
    )
}

missing_metadata <- setdiff(
    sample_names,
    samples_unique$sample
)

if (length(missing_metadata) > 0) {
    stop(
        "Samples found in counts but not in samplesheet: ",
        paste(missing_metadata, collapse = ", ")
    )
}


### ---------------------------
### 3. Read first featureCounts file
### ---------------------------

first_df <- read.delim(
    count_files[[1]],
    header = TRUE,
    comment.char = "#",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

colnames(first_df) <- trimws(colnames(first_df))
first_df$Geneid <- trimws(first_df$Geneid)

required_columns <- c("Geneid", "Length")

missing_columns <- setdiff(required_columns, colnames(first_df))

if (length(missing_columns) > 0) {
    stop(
        "Missing featureCounts columns in ",
        count_files[[1]],
        ": ",
        paste(missing_columns, collapse = ", ")
    )
}

gene_ids <- first_df$Geneid
gene_lengths <- first_df$Length
names(gene_lengths) <- gene_ids


### ---------------------------
### 4. Read counts
### ---------------------------

read_featurecounts <- function(file, expected_gene_ids) {

    df <- read.delim(
        file,
        header = TRUE,
        comment.char = "#",
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    colnames(df) <- trimws(colnames(df))
    df$Geneid <- trimws(df$Geneid)

    if (!"Geneid" %in% colnames(df)) {
        stop("Column 'Geneid' was not found in: ", file)
    }

    count_column <- colnames(df)[ncol(df)]

    values <- suppressWarnings(
        as.numeric(df[[count_column]])
    )

    if (anyNA(values)) {
        stop(
            "Non-numeric or missing counts found in: ",
            file,
            "; count column: ",
            count_column
        )
    }

    gene_positions <- match(expected_gene_ids, df$Geneid)

    if (anyNA(gene_positions)) {
        stop(
            "Gene IDs differ between featureCounts files. Problematic file: ",
            file
        )
    }

    values <- values[gene_positions]

    message(
        basename(file),
        ": ",
        sum(values),
        " assigned reads"
    )

    values
}

count_vectors <- lapply(
    count_files,
    read_featurecounts,
    expected_gene_ids = gene_ids
)

cts_fc <- do.call(cbind, count_vectors)

rownames(cts_fc) <- gene_ids
colnames(cts_fc) <- sample_names

storage.mode(cts_fc) <- "integer"


### ---------------------------
### 5. Validate matrix
### ---------------------------

count_totals <- colSums(cts_fc)

message("Count totals:")
print(count_totals)

zero_samples <- names(count_totals)[count_totals == 0]

if (length(zero_samples) > 0) {
    stop(
        "Samples with zero assigned counts: ",
        paste(zero_samples, collapse = ", "),
        ". Check the original featureCounts files and parameters."
    )
}

if (all(cts_fc == 0)) {
    stop("The complete featureCounts matrix contains only zeros.")
}

samples_unique <- samples_unique[
    colnames(cts_fc),
    ,
    drop = FALSE
]

if (!identical(colnames(cts_fc), rownames(samples_unique))) {
    stop("Count matrix columns and coldata row names are not aligned.")
}


### ---------------------------
### 6. Save result
### ---------------------------

saveRDS(
    list(
        counts = cts_fc,
        coldata = samples_unique,
        gene_lengths = gene_lengths
    ),
    snakemake@output[["raw_counts"]]
)