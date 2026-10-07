#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ============================================================
# Colours and GO categories
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

go_subcategory_colors <- c(
  "RNA processing" = "#4E79A7",
  "RNA splicing" = "#A0CBE8",
  "mRNA processing" = "#59A14F",
  "ribosome biogenesis" = "#8CD17D",
  "DNA damage" = "#F28E2B",
  "DNA repair" = "#FFBE7D",
  "response to DNA damage" = "#B07AA1",
  "response to stress" = "#E15759",
  "stress response" = "#FF9D9A",
  "oxidative stress" = "#9C755F",
  "heat shock" = "#BAB0AC",
  "chromatin" = "#76B7B2",
  "histone" = "#86BCB6",
  "nucleosome" = "#EDC948",
  "chromosome organization" = "#B6992D",
  "DNA methylation" = "#AF7AA1",
  "histone modification" = "#D4A6C8",
  "methyltransferase" = "#79706E",
  "acetyltransferase" = "#D7B5A6",
  "deacetylase" = "#9D7660",
  "epigenetic" = "#C7C7C7",
  "Other" = "grey80"
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

category_levels <- names(categories)

category_labels <- c(
  RNA_processing = "RNA\nprocessing",
  DNA_damage_response = "DNA damage\nresponse",
  Stress_response = "Stress\nresponse",
  Epigenetic_landscape = "Epigenetic\nlandscape"
)

# ============================================================
# Inputs and output directory
# ============================================================

outdir <- as.character(
  snakemake@output[["go_dir"]]
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

raw <- readRDS(
  snakemake@input[["raw_counts"]]
)

coldata <- raw$coldata
coldata$sample <- as.character(coldata$sample)
coldata$species <- as.character(coldata$species)
coldata$condition <- as.character(coldata$condition)

tpm <- as.matrix(
  readRDS(snakemake@input[["tpm"]])
)

storage.mode(tpm) <- "numeric"

remove_ids <- readLines(
  snakemake@input[["remove_ids"]],
  warn = FALSE
)

remove_ids <- unique(
  toupper(
    trimws(remove_ids)
  )
)

remove_ids <- remove_ids[
  !is.na(remove_ids) & nzchar(remove_ids)
]


# ============================================================
# Read the old-style final DEG annotation tables
# ============================================================

read_final_annotation <- function(path) {
  readr::read_delim(
    file = path,
    delim = ";",
    locale = readr::locale(
      decimal_mark = ","
    ),
    show_col_types = FALSE,
    progress = FALSE,
    trim_ws = TRUE,
    name_repair = "unique"
  ) %>%
    as.data.frame(
      stringsAsFactors = FALSE
    )
}

final_tables <- list(
  experimentalis = read_final_annotation(
    snakemake@input[["exp_final"]]
  ),
  gadabouti = read_final_annotation(
    snakemake@input[["gad_final"]]
  )
)

required_final_columns <- c(
  "Gene_Symbol",
  "log2FoldChange",
  "padj",
  "GO_Terms",
  "GO_Term_Description",
  "GTF_annotation",
  "Protein_ID"
)

prepare_final_annotation <- function(x, species_name) {
  missing_columns <- setdiff(
    required_final_columns,
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

  x <- x %>%
    mutate(
      GeneID = trimws(
        as.character(Gene_Symbol)
      ),
      gene_key = toupper(GeneID),
      gene_label = case_when(
        !is.na(GTF_annotation) &
          nzchar(trimws(GTF_annotation)) ~ paste0(
            GTF_annotation,
            " [",
            GeneID,
            "]"
          ),
        !is.na(Protein_ID) &
          nzchar(trimws(Protein_ID)) ~ paste0(
            Protein_ID,
            " [",
            GeneID,
            "]"
          ),
        TRUE ~ GeneID
      )
    )

  excluded <- x$gene_key %in% remove_ids

  message(
    species_name,
    ": ",
    nrow(x),
    " final DEGs loaded; ",
    sum(excluded),
    " excluded IDs removed; ",
    sum(!excluded),
    " retained"
  )

  x[
    !excluded,
    ,
    drop = FALSE
  ]
}

final_tables <- imap(
  final_tables,
  prepare_final_annotation
)

# ============================================================
# Identify GO-interest rows as in the previous notebook
# ============================================================

filter_go_interest <- function(df) {
  df <- as_tibble(df) %>%
    mutate(
      across(
        everything(),
        as.character
      )
    )

  purrr::imap_dfr(
    categories,
    function(patterns, category_name) {
      regex_pattern <- paste(
        patterns,
        collapse = "|"
      )

      df %>%
        filter(
          if_any(
            everything(),
            ~ str_detect(
              .x,
              regex(
                regex_pattern,
                ignore_case = TRUE
              )
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

assign_go_subcategory <- function(df) {
  if (!nrow(df)) {
    df$GO_subcategory <- character(0)
    return(df)
  }

  annotation_text <- paste(
    ifelse(
      is.na(df$GO_Term_Description),
      "",
      df$GO_Term_Description
    ),
    ifelse(
      is.na(df$GO_Terms),
      "",
      df$GO_Terms
    ),
    ifelse(
      is.na(df$GTF_annotation),
      "",
      df$GTF_annotation
    ),
    sep = " "
  )

  result <- rep(
    NA_character_,
    nrow(df)
  )

  for (category_name in names(categories)) {
    for (term in categories[[category_name]]) {
      matched <- str_detect(
        annotation_text,
        regex(
          term,
          ignore_case = TRUE
        )
      )

      result[
        matched & is.na(result)
      ] <- term
    }
  }

  result[is.na(result)] <- "Other"
  df$GO_subcategory <- result
  df
}

interest <- final_tables %>%
  map(filter_go_interest) %>%
  map(assign_go_subcategory)

for (species_name in names(interest)) {
  message(
    species_name,
    ": ",
    n_distinct(interest[[species_name]]$GeneID),
    " unique GO-interest DEGs"
  )

  readr::write_csv(
    interest[[species_name]],
    file.path(
      outdir,
      paste0(
        species_name,
        "_GO_interest.csv"
      )
    )
  )
}

message(
  "GO-interest DEG union: ",
  n_distinct(
    c(
      interest$experimentalis$GeneID,
      interest$gadabouti$GeneID
    )
  )
)

# ============================================================
# Heatmap metadata
# ============================================================

make_rows <- function(df) {
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

  df %>%
    transmute(
      gene_id = as.character(GeneID),
      label = as.character(gene_label),
      GO_category,
      GO_subcategory,
      DEG = "DEG"
    ) %>%
    distinct(
      gene_id,
      GO_category,
      .keep_all = TRUE
    )
}

rows <- map(
  interest,
  make_rows
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
      "No rows for ",
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

  row_metadata <- row_metadata %>%
    distinct(
      gene_id,
      .keep_all = TRUE
    )

  genes <- intersect(
    row_metadata$gene_id,
    rownames(tpm)
  )

  if (!length(samples) || !length(genes)) {
    message(
      "No overlapping TPM data for ",
      species_name
    )
    return(FALSE)
  }

  annotation_rows <- row_metadata[
    match(genes, row_metadata$gene_id),
    ,
    drop = FALSE
  ]

  expression_matrix <- tpm[
    genes,
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
    return(FALSE)
  }

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

    values <- as.numeric(plot_matrix)
    values <- values[is.finite(values)]

    breaks <- as.numeric(
      quantile(
        values,
        c(0, 0.5, 0.95),
        na.rm = TRUE
      )
    )

    if (length(unique(breaks)) < 3) {
      value_range <- range(
        values,
        na.rm = TRUE
      )
      if (
        !all(is.finite(value_range)) ||
        value_range[1] == value_range[2]
      ) {
        value_range <- c(0, 1)
      }
      breaks <- seq(
        value_range[1],
        value_range[2],
        length.out = 3
      )
    }

    color_function <- colorRamp2(
      breaks,
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

    plot_matrix[!is.finite(plot_matrix)] <- 0

    limit <- max(
      abs(plot_matrix),
      na.rm = TRUE
    )

    if (!is.finite(limit) || limit == 0) {
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

  condition_values <- coldata$condition[
    match(
      colnames(plot_matrix),
      coldata$sample
    )
  ]

  unknown_conditions <- setdiff(
    unique(condition_values),
    names(colors_condition)
  )

  if (length(unknown_conditions)) {
    stop(
      "Conditions without colours: ",
      paste(unknown_conditions, collapse = ", ")
    )
  }

  heatmap <- Heatmap(
    plot_matrix,
    name = legend_name,
    col = color_function,
    top_annotation = HeatmapAnnotation(
      Condition = condition_values,
      col = list(
        Condition = colors_condition
      )
    ),
    left_annotation = rowAnnotation(
      GO = annotation_rows$GO_category,
      DEG = annotation_rows$DEG,
      col = list(
        GO = go_colors,
        DEG = deg_colors
      )
    ),
    row_labels = annotation_rows$label,
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 10),
    column_title = paste(
      "GO-interest DEGs",
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
# Generate heatmaps
# ============================================================

for (species_name in names(rows)) {
  species_rows <- rows[[species_name]]

  for (mode in c("logTPM", "zscore_logTPM")) {
    plot_heat(
      species_name,
      species_rows,
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

  for (category_name in category_levels) {
    category_rows <- species_rows[
      species_rows$GO_category == category_name,
      ,
      drop = FALSE
    ]

    for (mode in c("logTPM", "zscore_logTPM")) {
      plot_heat(
        species_name,
        category_rows,
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

# ============================================================
# GO-interest summary and barplots
# ============================================================

summary_table <- imap_dfr(
  interest,
  function(df, species_name) {
    df %>%
      distinct(
        GeneID,
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
  }
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

make_summary_plot <- function(plot_data, title_text) {
  ggplot(
    plot_data,
    aes(
      x = factor(
        GO_category,
        levels = category_levels
      ),
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
      size = 3,
      color = "black"
    ) +
    scale_fill_manual(
      values = go_subcategory_colors,
      na.value = "grey80",
      drop = FALSE
    ) +
    scale_x_discrete(
      labels = category_labels
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = title_text,
      x = "GO category",
      y = "Number of DE genes",
      fill = "GO term group"
    )
}

combined_plot <- make_summary_plot(
  summary_table,
  "GO-interest categories among DEGs"
) +
  facet_wrap(
    ~ species
  )

ggsave(
  file.path(
    summary_dir,
    "GO_interest_DEG_barplot.pdf"
  ),
  plot = combined_plot,
  width = 12,
  height = 6,
  dpi = 300
)

for (species_name in unique(summary_table$species)) {
  species_summary <- summary_table %>%
    filter(
      species == species_name
    )

  species_plot <- make_summary_plot(
    species_summary,
    paste(
      "GO-interest categories —",
      species_name
    )
  )

  ggsave(
    file.path(
      summary_dir,
      paste0(
        "GO_interest_DEG_barplot_",
        species_name,
        ".pdf"
      )
    ),
    plot = species_plot,
    width = 10,
    height = 6,
    dpi = 300
  )
}

message(
  "GO-interest analysis completed successfully."
)