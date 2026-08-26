#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
})

colors_condition <- c(
  control = "#1b9e77",
  anhydrobiosis = "#d95f02"
)

to_matrix <- function(x) {
  if (is.matrix(x)) {
    return(x)
  }

  as.matrix(
    SummarizedExperiment::assay(x)
  )
}

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
    "coldata must contain a 'sample' column. Found: ",
    paste(colnames(coldata), collapse = ", ")
  )
}

centered <- list()
all_values <- numeric()

for (sp in names(res_list)) {
  if (!sp %in% names(vst_list)) {
    stop(
      "Species '",
      sp,
      "' was not found in vst_list. Available species: ",
      paste(names(vst_list), collapse = ", ")
    )
  }

  mat <- to_matrix(
    vst_list[[sp]]
  )

  species_results <- res_list[[sp]]

  sig_genes <- species_results %>%
    filter(
      !is.na(sig),
      sig == "yes"
    ) %>%
    pull(GeneID)

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

  # Remove genes with zero variance before calculating Z-scores
  gene_sd <- apply(
    expression_matrix,
    1,
    sd,
    na.rm = TRUE
  )

  keep <- is.finite(gene_sd) & gene_sd > 0

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

  all_values <- c(
    all_values,
    as.numeric(z)
  )

  message(
    sp,
    ": ",
    nrow(z),
    " significant genes included in the heatmap"
  )
}

max_abs <- max(
  abs(all_values),
  na.rm = TRUE
)

if (!is.finite(max_abs) || max_abs == 0) {
  max_abs <- 1
}

col_fun <- colorRamp2(
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

outputs <- setNames(
  as.character(snakemake@output),
  c(
    "experimentalis",
    "gadabouti"
  )
)

for (sp in names(outputs)) {
  if (!sp %in% names(centered)) {
    stop(
      "No heatmap matrix was generated for ",
      sp
    )
  }

  z <- centered[[sp]]

  sample_match <- match(
    colnames(z),
    coldata$sample
  )

  if (anyNA(sample_match)) {
    stop(
      "Samples from the VST matrix were not found in coldata: ",
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

  heatmap <- Heatmap(
    z,
    name = "Z-score",
    col = col_fun,

    top_annotation = HeatmapAnnotation(
      df = annotation_data,
      col = list(
        condition = colors_condition
      ),
      annotation_legend_param = list(
        condition = list(
          title = "Condition"
        )
      )
    ),

    show_row_names = FALSE,
    show_column_names = TRUE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,

    column_title = paste0(
      "Heatmap ",
      sp,
      " — ",
      nrow(z),
      " significant DEGs"
    ),

    column_names_gp = gpar(
      fontsize = 12
    ),

    heatmap_legend_param = list(
      title = "Z-score",
      labels_gp = gpar(
        fontsize = 10
      ),
      title_gp = gpar(
        fontsize = 12,
        fontface = "bold"
      )
    )
  )

  pdf(
    outputs[[sp]],
    width = 8,
    height = 8
  )

  draw(heatmap)
  dev.off()
}