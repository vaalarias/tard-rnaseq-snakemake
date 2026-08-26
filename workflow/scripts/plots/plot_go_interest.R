#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ============================================================
# Configuration
# ============================================================

colors_condition <- c(
  control = "#1b9e77",
  anhydrobiosis = "#d95f02",
  rehydrated = "#d95f02"
)

go_colors <- c(
  RNA_processing = "#1b9e77",
  DNA_damage_response = "#7570b3",
  Stress_response = "#e7298a",
  Epigenetic_landscape = "#66a61e",
  Other = "gray85"
)

deg_colors <- c(
  DEG = "#d73027",
  not_DEG = "gray85"
)

categories <- list(
  RNA_processing = c(
    "RNA processing",
    "RNA splicing",
    "mRNA processing",
    "ribosome biogenesis"
  ),

  DNA_damage_response = c(
    "DNA damage",
    "DNA repair",
    "response to DNA damage"
  ),

  Stress_response = c(
    "response to stress",
    "stress response",
    "oxidative stress",
    "heat shock"
  ),

  Epigenetic_landscape = c(
    "chromatin",
    "histone",
    "nucleosome",
    "chromosome organization",
    "DNA methylation",
    "histone modification",
    "methyltransferase",
    "acetyltransferase",
    "deacetylase",
    "epigenetic"
  )
)

# ============================================================
# Read inputs
# ============================================================

outdir <- snakemake@output[["go_dir"]]

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

raw <- readRDS(
  snakemake@input[["raw_counts"]]
)

coldata <- raw$coldata

tpm <- as.matrix(
  readRDS(snakemake@input[["tpm"]])
)

storage.mode(tpm) <- "numeric"

res <- readRDS(
  snakemake@input[["deseq_results"]]
)

annot_exp <- readRDS(
  snakemake@input[["annot_exp"]]
)

annot_gad <- readRDS(
  snakemake@input[["annot_gad"]]
)

term2gene <- as.data.frame(
  readRDS(snakemake@input[["term2gene"]]),
  stringsAsFactors = FALSE
)

term2name <- as.data.frame(
  readRDS(snakemake@input[["term2name"]]),
  stringsAsFactors = FALSE
)

required_coldata <- c(
  "sample",
  "species",
  "condition"
)

missing_coldata <- setdiff(
  required_coldata,
  colnames(coldata)
)

if (length(missing_coldata)) {
  stop(
    "Missing coldata columns: ",
    paste(missing_coldata, collapse = ", ")
  )
}

if (!"GeneID" %in% colnames(annot_exp)) {
  stop(
    "experimentalis_annotated.rds has no GeneID column."
  )
}

if (!"GeneID" %in% colnames(annot_gad)) {
  stop(
    "gadabouti_annotated.rds has no GeneID column."
  )
}

if (ncol(term2gene) < 2) {
  stop(
    "TERM2GENE.rds must contain at least two columns."
  )
}

if (ncol(term2name) < 2) {
  stop(
    "TERM2NAME.rds must contain at least two columns."
  )
}

# ============================================================
# General helpers
# ============================================================

clean_gene_key <- function(x) {
  x <- trimws(
    as.character(x)
  )

  x <- sub(
    "^LOC",
    "",
    x,
    ignore.case = TRUE
  )

  toupper(x)
}

# ============================================================
# Read GAF and recover gene-product labels
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
  snakemake@input[["gaf"]],
  comment.char = "!",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  fill = TRUE,
  check.names = FALSE
)

if (ncol(gaf) < 15) {
  stop(
    "Invalid GAF: expected at least 15 columns, found ",
    ncol(gaf)
  )
}

if (ncol(gaf) > length(gaf_columns)) {
  stop(
    "Invalid GAF: expected no more than ",
    length(gaf_columns),
    " columns, found ",
    ncol(gaf)
  )
}

colnames(gaf) <- gaf_columns[
  seq_len(ncol(gaf))
]

gaf <- gaf %>%
  mutate(
    across(
      c(
        DB_Object_ID,
        DB_Object_Symbol,
        DB_Object_Name,
        Taxon
      ),
      as.character
    ),

    gene_product = case_when(
      !is.na(DB_Object_Name) &
        nzchar(trimws(DB_Object_Name)) &
        DB_Object_Name != "-" ~ DB_Object_Name,

      !is.na(DB_Object_Symbol) &
        nzchar(trimws(DB_Object_Symbol)) &
        DB_Object_Symbol != "-" ~ DB_Object_Symbol,

      TRUE ~ as.character(DB_Object_ID)
    )
  )
  
gaf_gene_labels <- bind_rows(
  gaf %>%
    transmute(
      gene_key = clean_gene_key(DB_Object_ID),
      gene_product
    ),

  gaf %>%
    transmute(
      gene_key = clean_gene_key(DB_Object_Symbol),
      gene_product
    )
) %>%
  filter(
    !is.na(gene_key),
    nzchar(gene_key),
    !is.na(gene_product),
    nzchar(gene_product),
    gene_product != "-"
  ) %>%
  distinct(
    gene_key,
    .keep_all = TRUE
  )

message(
  "Gene labels recovered from GAF: ",
  nrow(gaf_gene_labels)
)

# ============================================================
# Save GAF taxon summary
# ============================================================

taxon <- sub(
  "taxon:",
  "",
  gaf$Taxon
)

taxon_summary <- as.data.frame(
  table(taxon),
  stringsAsFactors = FALSE
)

colnames(taxon_summary) <- c(
  "Taxon",
  "n"
)

taxon_summary$percent <- (
  taxon_summary$n /
    sum(taxon_summary$n)
) * 100

taxon_summary <- taxon_summary[
  order(-taxon_summary$n),
  ,
  drop = FALSE
]

write.csv(
  taxon_summary,
  file.path(
    outdir,
    "GAF_taxon_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# Prepare TERM2GENE and TERM2NAME
# ============================================================

term2gene <- term2gene[
  ,
  1:2,
  drop = FALSE
]

colnames(term2gene) <- c(
  "GO_ID",
  "annotation_gene_id"
)

term2name <- term2name[
  ,
  1:2,
  drop = FALSE
]

colnames(term2name) <- c(
  "GO_ID",
  "GO_Term_Description"
)

term2gene <- term2gene %>%
  mutate(
    GO_ID = as.character(GO_ID),

    annotation_gene_id = as.character(
      annotation_gene_id
    ),

    gene_key = clean_gene_key(
      annotation_gene_id
    )
  ) %>%
  filter(
    !is.na(GO_ID),
    nzchar(GO_ID),
    !is.na(gene_key),
    nzchar(gene_key)
  ) %>%
  distinct(
    GO_ID,
    gene_key,
    .keep_all = TRUE
  )

term2name <- term2name %>%
  mutate(
    GO_ID = as.character(GO_ID),

    GO_Term_Description = as.character(
      GO_Term_Description
    )
  ) %>%
  filter(
    !is.na(GO_ID),
    nzchar(GO_ID)
  ) %>%
  distinct(
    GO_ID,
    .keep_all = TRUE
  )

go_annotations <- term2gene %>%
  left_join(
    term2name,
    by = "GO_ID"
  )

# ============================================================
# Match annotation GeneIDs against GO mappings
# ============================================================

make_gene_lookup <- function(annotation_table) {
  annotation_table %>%
    transmute(
      GeneID = as.character(GeneID),

      gene_key = clean_gene_key(
        GeneID
      )
    ) %>%
    filter(
      !is.na(GeneID),
      nzchar(GeneID),
      !is.na(gene_key),
      nzchar(gene_key)
    ) %>%
    distinct(
      gene_key,
      .keep_all = TRUE
    )
}

lookup_exp <- make_gene_lookup(
  annot_exp
)

lookup_gad <- make_gene_lookup(
  annot_gad
)

make_species_annotation <- function(gene_lookup) {
  go_annotations %>%
    inner_join(
      gene_lookup,
      by = "gene_key"
    ) %>%
    left_join(
      gaf_gene_labels,
      by = "gene_key"
    ) %>%
    mutate(
      gene_label = case_when(
        !is.na(gene_product) &
          nzchar(gene_product) ~ paste0(
            gene_product,
            " [",
            GeneID,
            "]"
          ),

        TRUE ~ GeneID
      )
    ) %>%
    select(
      GeneID,
      GO_ID,
      GO_Term_Description,
      annotation_gene_id,
      gene_product,
      gene_label
    ) %>%
    distinct()
}

go_exp <- make_species_annotation(
  lookup_exp
)

go_gad <- make_species_annotation(
  lookup_gad
)

message(
  "GO-annotated experimentalis genes: ",
  n_distinct(go_exp$GeneID)
)

message(
  "GO-annotated gadabouti genes: ",
  n_distinct(go_gad$GeneID)
)

# ============================================================
# Select GO categories of interest
# ============================================================

filter_interest <- function(df) {
  if (!nrow(df)) {
    df$GO_category <- character(0)
    return(df)
  }

  imap_dfr(
    categories,
    function(patterns, category_name) {
      category_pattern <- paste(
        patterns,
        collapse = "|"
      )

      df %>%
        filter(
          !is.na(GO_Term_Description),

          str_detect(
            GO_Term_Description,
            regex(
              category_pattern,
              ignore_case = TRUE
            )
          )
        ) %>%
        mutate(
          GO_category = category_name
        )
    }
  ) %>%
    distinct()
}

interest <- list(
  experimentalis = filter_interest(
    go_exp
  ),

  gadabouti = filter_interest(
    go_gad
  )
)

# ============================================================
# Assign GO subcategories
# ============================================================

assign_subcategory <- function(df) {
  if (!nrow(df)) {
    df$GO_subcategory <- character(0)
    return(df)
  }

  description <- as.character(
    df$GO_Term_Description
  )

  result <- rep(
    "Other",
    nrow(df)
  )

  for (category_name in names(categories)) {
    for (term in categories[[category_name]]) {
      matched <- str_detect(
        description,
        regex(
          term,
          ignore_case = TRUE
        )
      )

      result[
        matched &
          result == "Other"
      ] <- term
    }
  }

  df$GO_subcategory <- result
  df
}

interest <- lapply(
  interest,
  assign_subcategory
)

# ============================================================
# Save GO-interest tables with gene labels
# ============================================================

write_csv(
  interest$experimentalis,
  file.path(
    outdir,
    "experimentalis_GO_interest.csv"
  )
)

write_csv(
  interest$gadabouti,
  file.path(
    outdir,
    "gadabouti_GO_interest.csv"
  )
)

# ============================================================
# Construct heatmap row metadata
# ============================================================

make_rows <- function(df, species) {
  if (!nrow(df)) {
    return(
      tibble(
        gene_id = character(),
        label = character(),
        GO_category = character(),
        GO_subcategory = character(),
        DEG = character()
      )
    )
  }

  x <- df %>%
    transmute(
      gene_id = as.character(GeneID),

      label = case_when(
        !is.na(gene_label) &
          nzchar(gene_label) ~ gene_label,

        TRUE ~ as.character(GeneID)
      ),

      GO_category,
      GO_subcategory
    ) %>%
    distinct(
      gene_id,
      GO_category,
      .keep_all = TRUE
    )

  species_results <- res[[species]]

  significant_genes <- species_results$GeneID[
    !is.na(species_results$sig) &
      species_results$sig == "yes"
  ]

  x$DEG <- ifelse(
    x$gene_id %in% significant_genes,
    "DEG",
    "not_DEG"
  )

  x
}

rows <- list(
  experimentalis = make_rows(
    interest$experimentalis,
    "experimentalis"
  ),

  gadabouti = make_rows(
    interest$gadabouti,
    "gadabouti"
  )
)

# ============================================================
# Heatmap function
# ============================================================

plot_heat <- function(
  species_name,
  row_metadata,
  mode,
  output_file
) {
  if (!nrow(row_metadata)) {
    message(
      "No GO-interest rows for ",
      species_name,
      " | ",
      mode
    )

    return(FALSE)
  }

  samples <- coldata$sample[
    coldata$species == species_name
  ]

  samples <- intersect(
    samples,
    colnames(tpm)
  )

  if (!length(samples)) {
    message(
      "No TPM samples for ",
      species_name
    )

    return(FALSE)
  }

  # One row per gene in each heatmap.
  row_metadata <- row_metadata %>%
    distinct(
      gene_id,
      .keep_all = TRUE
    )

  genes <- intersect(
    row_metadata$gene_id,
    rownames(tpm)
  )

  if (!length(genes)) {
    message(
      "No GO-interest genes overlap TPM for ",
      species_name
    )

    return(FALSE)
  }

  annotation_rows <- row_metadata[
    match(
      genes,
      row_metadata$gene_id
    ),
    ,
    drop = FALSE
  ]

  expression_matrix <- tpm[
    annotation_rows$gene_id,
    samples,
    drop = FALSE
  ]

  keep <- rowSums(
    expression_matrix > 0,
    na.rm = TRUE
  ) > 0

  expression_matrix <- expression_matrix[
    keep,
    ,
    drop = FALSE
  ]

  annotation_rows <- annotation_rows[
    keep,
    ,
    drop = FALSE
  ]

  if (!nrow(expression_matrix)) {
    message(
      "All GO-interest genes have zero TPM for ",
      species_name
    )

    return(FALSE)
  }

  category_levels <- c(
    names(categories),
    "Other"
  )

  annotation_rows$GO_category <- factor(
    annotation_rows$GO_category,
    levels = category_levels
  )

  ordering <- order(
    annotation_rows$GO_category,
    annotation_rows$label,
    annotation_rows$gene_id
  )

  expression_matrix <- expression_matrix[
    ordering,
    ,
    drop = FALSE
  ]

  annotation_rows <- annotation_rows[
    ordering,
    ,
    drop = FALSE
  ]

  if (mode == "logTPM") {
    plot_matrix <- log10(
      expression_matrix + 1
    )

    values <- as.numeric(
      plot_matrix
    )

    values <- values[
      is.finite(values)
    ]

    quantiles <- quantile(
      values,
      c(0, 0.5, 0.95),
      na.rm = TRUE
    )

    if (length(unique(quantiles)) < 3) {
      value_range <- range(
        values,
        na.rm = TRUE
      )

      if (
        !all(is.finite(value_range)) ||
          value_range[1] == value_range[2]
      ) {
        value_range <- c(
          0,
          1
        )
      }

      quantiles <- seq(
        value_range[1],
        value_range[2],
        length.out = 3
      )
    }

    color_function <- colorRamp2(
      quantiles,
      c(
        "white",
        "#fee08b",
        "#d73027"
      )
    )

    legend_name <- "log10(TPM+1)"
  } else if (mode == "zscore_logTPM") {
    plot_matrix <- t(
      scale(
        t(
          log10(expression_matrix + 1)
        )
      )
    )

    plot_matrix[
      !is.finite(plot_matrix)
    ] <- 0

    limit <- max(
      abs(plot_matrix),
      na.rm = TRUE
    )

    if (
      !is.finite(limit) ||
        limit == 0
    ) {
      limit <- 1
    }

    color_function <- colorRamp2(
      c(-limit, 0, limit),
      c(
        "#2166ac",
        "white",
        "#b2182b"
      )
    )

    legend_name <- "Z-score"
  } else {
    stop(
      "Unsupported heatmap mode: ",
      mode
    )
  }

  condition_values <- as.character(
    coldata$condition[
      match(
        colnames(plot_matrix),
        coldata$sample
      )
    ]
  )

  unknown_conditions <- setdiff(
    unique(condition_values),
    names(colors_condition)
  )

  if (length(unknown_conditions)) {
    stop(
      "Conditions without colors: ",
      paste(
        unknown_conditions,
        collapse = ", "
      )
    )
  }

  top_annotation <- HeatmapAnnotation(
    Condition = condition_values,
    col = list(
      Condition = colors_condition
    )
  )

  left_annotation <- rowAnnotation(
    GO = annotation_rows$GO_category,
    DEG = annotation_rows$DEG,
    col = list(
      GO = go_colors,
      DEG = deg_colors
    )
  )

  heatmap <- Heatmap(
    plot_matrix,
    name = legend_name,
    col = color_function,
    top_annotation = top_annotation,
    left_annotation = left_annotation,

    # Gene product/annotation, not GO description.
    row_labels = annotation_rows$label,

    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_rows = FALSE,
    cluster_columns = TRUE,

    row_names_gp = gpar(
      fontsize = 6
    ),

    column_title = paste(
      "GO interest genes",
      species_name,
      mode,
      sep = " | "
    )
  )

  dir.create(
    dirname(output_file),
    recursive = TRUE,
    showWarnings = FALSE
  )

  pdf(
    output_file,
    width = 8,
    height = max(
      5,
      min(
        14,
        nrow(plot_matrix) * 0.18
      )
    )
  )

  draw(heatmap)
  dev.off()

  message(
    "Saved: ",
    output_file
  )

  TRUE
}

# ============================================================
# Generate heatmaps and summaries
# ============================================================

summary_list <- list()

for (species_name in names(rows)) {
  deg_rows <- rows[[species_name]][
    rows[[species_name]]$DEG == "DEG",
    ,
    drop = FALSE
  ]

  summary_list[[species_name]] <- deg_rows %>%
    distinct(
      gene_id,
      GO_category,
      GO_subcategory
    ) %>%
    count(
      GO_category,
      GO_subcategory,
      name = "n"
    ) %>%
    mutate(
      species = species_name
    )

  for (category_name in names(categories)) {
    category_rows <- rows[[species_name]][
      rows[[species_name]]$GO_category ==
        category_name,
      ,
      drop = FALSE
    ]

    for (mode in c(
      "logTPM",
      "zscore_logTPM"
    )) {
      plot_heat(
        species_name,
        category_rows,
        mode,
        file.path(
          outdir,
          paste0(
            species_name,
            "_",
            category_name,
            "_",
            mode,
            ".pdf"
          )
        )
      )

      category_deg_rows <- category_rows[
        category_rows$DEG == "DEG",
        ,
        drop = FALSE
      ]

      if (nrow(category_deg_rows)) {
        plot_heat(
          species_name,
          category_deg_rows,
          mode,
          file.path(
            outdir,
            "DEG_only",
            "by_category",
            category_name,
            paste0(
              species_name,
              "_",
              mode,
              ".pdf"
            )
          )
        )
      }
    }
  }

  for (mode in c(
    "logTPM",
    "zscore_logTPM"
  )) {
    if (nrow(deg_rows)) {
      plot_heat(
        species_name,
        deg_rows,
        mode,
        file.path(
          outdir,
          "DEG_only",
          paste0(
            species_name,
            "_GO_interest_DEG_only_",
            mode,
            ".pdf"
          )
        )
      )
    }
  }
}

# ============================================================
# GO DEG summary
# ============================================================

summary_table <- bind_rows(
  summary_list
)

summary_dir <- file.path(
  outdir,
  "GO_summary"
)

dir.create(
  summary_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  summary_table,
  file.path(
    summary_dir,
    "GO_interest_DEG_summary.csv"
  ),
  row.names = FALSE
)

if (nrow(summary_table)) {
  summary_plot <- ggplot(
    summary_table,
    aes(
      x = GO_category,
      y = n,
      fill = GO_subcategory
    )
  ) +
    geom_col(
      color = "grey25",
      linewidth = 0.2
    ) +
    geom_text(
      aes(label = n),
      position = position_stack(
        vjust = 0.5
      ),
      size = 3
    ) +
    facet_wrap(
      ~ species
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      axis.text.x = element_text(
        angle = 35,
        hjust = 1
      ),

      panel.grid.major.x = element_blank()
    ) +
    labs(
      title = "GO interest categories among DEGs",
      x = "GO category",
      y = "Number of DEG genes",
      fill = "GO term group"
    )

  ggsave(
    file.path(
      summary_dir,
      "GO_interest_DEG_barplot.pdf"
    ),
    summary_plot,
    width = 12,
    height = 6,
    dpi = 300
  )
} else {
  pdf(
    file.path(
      summary_dir,
      "GO_interest_DEG_barplot.pdf"
    ),
    width = 8,
    height = 5
  )

  plot.new()

  title(
    main = "No GO interest categories found among DEGs"
  )

  dev.off()
}