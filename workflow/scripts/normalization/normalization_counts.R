#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
})

# -------------------------------------------------------------------------
# Load and validate input
# -------------------------------------------------------------------------

raw <- readRDS(
  snakemake@input[["raw_counts"]]
)

required_objects <- c(
  "counts",
  "coldata",
  "gene_lengths"
)

missing_objects <- setdiff(
  required_objects,
  names(raw)
)

if (length(missing_objects) > 0) {
  stop(
    "Missing objects in raw-count RDS: ",
    paste(missing_objects, collapse = ", ")
  )
}

cts <- as.matrix(raw$counts)
coldata <- as.data.frame(raw$coldata)
gene_lengths <- raw$gene_lengths

if (is.null(rownames(cts))) {
  stop("Count matrix has no gene row names.")
}

if (is.null(colnames(cts))) {
  stop("Count matrix has no sample column names.")
}

if (is.null(rownames(coldata))) {
  stop("coldata must have sample IDs as row names.")
}

if (!setequal(colnames(cts), rownames(coldata))) {
  stop(
    "Samples in count matrix and coldata do not match."
  )
}

# Put metadata in exactly the same order as count columns.
coldata <- coldata[
  colnames(cts),
  ,
  drop = FALSE
]

if (
  anyNA(cts) ||
  any(!is.finite(cts)) ||
  any(cts < 0)
) {
  stop(
    "Counts contain NA, non-finite, or negative values."
  )
}

if (any(abs(cts - round(cts)) > 1e-8)) {
  stop(
    "DESeq2 requires integer counts, but non-integer ",
    "values were found."
  )
}

cts <- round(cts)
storage.mode(cts) <- "integer"

# -------------------------------------------------------------------------
# Match gene lengths to the count matrix
# -------------------------------------------------------------------------

if (
  is.matrix(gene_lengths) ||
  is.data.frame(gene_lengths)
) {
  if (ncol(gene_lengths) != 1) {
    stop(
      "gene_lengths must be a vector or a one-column object."
    )
  }

  gene_lengths <- gene_lengths[, 1]
}

if (!is.null(names(gene_lengths))) {
  missing_lengths <- setdiff(
    rownames(cts),
    names(gene_lengths)
  )

  if (length(missing_lengths) > 0) {
    stop(
      "Gene lengths are missing for ",
      length(missing_lengths),
      " genes. Examples: ",
      paste(head(missing_lengths), collapse = ", ")
    )
  }

  gene_lengths <- gene_lengths[
    rownames(cts)
  ]

} else if (length(gene_lengths) != nrow(cts)) {
  stop(
    "Unnamed gene_lengths does not have the same ",
    "length as the number of genes."
  )
}

gene_lengths <- as.numeric(
  gene_lengths
)

names(gene_lengths) <- rownames(cts)

if (
  anyNA(gene_lengths) ||
  any(!is.finite(gene_lengths)) ||
  any(gene_lengths <= 0)
) {
  stop(
    "Gene lengths contain NA, non-finite, zero, ",
    "or negative values."
  )
}

# -------------------------------------------------------------------------
# DESeq2 object and size factors
# -------------------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design = ~1
)

dds <- estimateSizeFactors(
  dds
)

# -------------------------------------------------------------------------
# VST matrix
# -------------------------------------------------------------------------

vst_object <- vst(
  dds,
  blind = TRUE
)

vst_matrix <- assay(
  vst_object
)

saveRDS(
  vst_matrix,
  snakemake@output[["vst"]]
)

# -------------------------------------------------------------------------
# rlog matrix
# -------------------------------------------------------------------------

rlog_object <- rlog(
  dds,
  blind = TRUE
)

rlog_matrix <- assay(
  rlog_object
)

saveRDS(
  rlog_matrix,
  snakemake@output[["rlog"]]
)

# -------------------------------------------------------------------------
# TMM-normalized CPM
# -------------------------------------------------------------------------

dge <- DGEList(
  counts = cts
)

dge <- calcNormFactors(
  dge,
  method = "TMM"
)

cpm_matrix <- cpm(
  dge,
  log = FALSE,
  normalized.lib.sizes = TRUE
)

saveRDS(
  cpm_matrix,
  snakemake@output[["cpm"]]
)

# -------------------------------------------------------------------------
# TPM
# -------------------------------------------------------------------------

gene_lengths_kb <- (
  gene_lengths /
  1000
)

rpk_matrix <- sweep(
  cts,
  MARGIN = 1,
  STATS = gene_lengths_kb,
  FUN = "/"
)

rpk_sums <- colSums(
  rpk_matrix
)

if (any(rpk_sums <= 0)) {
  stop(
    "At least one sample has a non-positive total RPK."
  )
}

tpm_matrix <- sweep(
  rpk_matrix,
  MARGIN = 2,
  STATS = rpk_sums,
  FUN = "/"
) * 1e6

tpm_column_sums <- colSums(
  tpm_matrix
)

if (
  !all(
    abs(tpm_column_sums - 1e6) < 1
  )
) {
  stop(
    "TPM validation failed. Column sums: ",
    paste(
      round(tpm_column_sums, 2),
      collapse = ", "
    )
  )
}

saveRDS(
  tpm_matrix,
  snakemake@output[["tpm"]]
)