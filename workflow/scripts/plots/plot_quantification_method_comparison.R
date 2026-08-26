#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(tximport)
  library(rtracklayer)
  library(patchwork)
  library(scales)
})

# -------------------------------------------------------------------------
# Inputs and outputs
# -------------------------------------------------------------------------

samplesheet_path <- snakemake@input[["samplesheet"]]
gtf_path <- snakemake@input[["gtf"]]

salmon_files <- unlist(
  snakemake@input[["salmon"]]
)

featurecounts_files <- unlist(
  snakemake@input[["featurecounts"]]
)

sample_ids <- unlist(
  snakemake@params[["samples"]]
)

method_plot_output <- snakemake@output[["method_plot"]]
overlap_plot_output <- snakemake@output[["overlap_plot"]]
combined_plot_output <- snakemake@output[["combined_plot"]]
sample_table_output <- snakemake@output[["sample_table"]]
overlap_table_output <- snakemake@output[["overlap_table"]]

dir.create(
  dirname(combined_plot_output),
  recursive = TRUE,
  showWarnings = FALSE
)

if (
  length(sample_ids) != length(salmon_files) ||
  length(sample_ids) != length(featurecounts_files)
) {
  stop(
    "Sample IDs, Salmon files, and featureCounts files ",
    "must have the same length."
  )
}

names(salmon_files) <- sample_ids
names(featurecounts_files) <- sample_ids

missing_salmon <- salmon_files[!file.exists(salmon_files)]
missing_featurecounts <- featurecounts_files[
  !file.exists(featurecounts_files)
]

if (length(missing_salmon) > 0) {
  stop(
    "Missing Salmon files:\n",
    paste(missing_salmon, collapse = "\n")
  )
}

if (length(missing_featurecounts) > 0) {
  stop(
    "Missing featureCounts files:\n",
    paste(missing_featurecounts, collapse = "\n")
  )
}

# -------------------------------------------------------------------------
# Metadata
# -------------------------------------------------------------------------

metadata <- read_csv(
  samplesheet_path,
  comment = "#",
  show_col_types = FALSE
) %>%
  filter(sample %in% sample_ids) %>%
  distinct(sample, .keep_all = TRUE) %>%
  mutate(
    species = factor(
      species,
      levels = c(
        "experimentalis",
        "gadabouti"
      ),
      labels = c(
        "Pam. experimentalis",
        "Pam. gadabouti"
      )
    ),
    condition = factor(
      condition,
      levels = c(
        "control",
        "anhydrobiosis"
      ),
      labels = c(
        "Control",
        "Post-anhydrobiosis"
      )
    )
  )

if (!setequal(metadata$sample, sample_ids)) {
  stop(
    "The samplesheet does not contain exactly the requested samples."
  )
}

metadata <- metadata %>%
  arrange(match(sample, sample_ids))

# -------------------------------------------------------------------------
# Transcript-to-gene mapping
# -------------------------------------------------------------------------

gtf <- import(gtf_path)

tx2gene <- as.data.frame(gtf) %>%
  filter(
    type == "transcript",
    !is.na(transcript_id),
    !is.na(gene_id)
  ) %>%
  transmute(
    TXNAME = as.character(transcript_id),
    GENEID = as.character(gene_id)
  ) %>%
  distinct()

if (nrow(tx2gene) == 0) {
  stop(
    "No transcript_id-to-gene_id relationships were found in the GTF."
  )
}

# Check compatibility between the Salmon index and GTF.
first_quant <- read_tsv(
  salmon_files[[1]],
  show_col_types = FALSE
)

mapping_fraction <- mean(
  first_quant$Name %in% tx2gene$TXNAME
)

message(
  "Fraction of Salmon transcript IDs found in GTF: ",
  round(mapping_fraction, 4)
)

if (mapping_fraction < 0.90) {
  stop(
    "Less than 90% of Salmon transcript IDs match transcript_id ",
    "values in the GTF. Confirm that Salmon and featureCounts used ",
    "the same annotation."
  )
}

# -------------------------------------------------------------------------
# Salmon: transcript-level estimates aggregated to genes
# -------------------------------------------------------------------------

salmon_import <- tximport(
  files = salmon_files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "no",
  ignoreTxVersion = FALSE
)

salmon_detection <- as.data.frame(
  salmon_import$counts
) %>%
  tibble::rownames_to_column("gene_id") %>%
  pivot_longer(
    cols = -gene_id,
    names_to = "sample",
    values_to = "count"
  ) %>%
  transmute(
    sample,
    gene_id,
    method = "Salmon",
    count,
    detected = count > 0
  )

# -------------------------------------------------------------------------
# STAR-featureCounts
# -------------------------------------------------------------------------

read_featurecounts <- function(path, sample) {
  counts <- read_tsv(
    path,
    comment = "#",
    show_col_types = FALSE,
    name_repair = "minimal"
  )

  if (!"Geneid" %in% colnames(counts)) {
    stop(
      "Geneid column not found in featureCounts file: ",
      path
    )
  }

  count_column <- tail(colnames(counts), 1)

  counts %>%
    transmute(
      sample = sample,
      gene_id = as.character(Geneid),
      method = "STAR-featureCounts",
      count = as.numeric(.data[[count_column]]),
      detected = count > 0
    )
}

featurecounts_detection <- bind_rows(
  Map(
    read_featurecounts,
    featurecounts_files,
    names(featurecounts_files)
  )
)

# -------------------------------------------------------------------------
# Per-sample detection
# -------------------------------------------------------------------------

detection <- bind_rows(
  featurecounts_detection,
  salmon_detection
) %>%
  inner_join(
    metadata %>%
      select(
        sample,
        species,
        condition,
        treatment
      ),
    by = "sample"
  )

sample_summary <- detection %>%
  group_by(
    sample,
    species,
    condition,
    treatment,
    method
  ) %>%
  summarise(
    detected_features = sum(detected),
    .groups = "drop"
  )

write_tsv(
  sample_summary,
  sample_table_output
)

plot_summary <- sample_summary %>%
  group_by(
    species,
    condition,
    method
  ) %>%
  summarise(
    n_samples = n(),
    mean_detected = mean(detected_features),
    sd_detected = sd(detected_features),
    .groups = "drop"
  )

method_colors <- c(
  "STAR-featureCounts" = "#4B2588",
  "Salmon" = "#D98C3F"
)

method_plot <- ggplot(
  plot_summary,
  aes(
    x = condition,
    y = mean_detected,
    fill = method
  )
) +
  geom_col(
    position = position_dodge(width = 0.78),
    width = 0.68,
    colour = "black",
    linewidth = 0.25
  ) +
  geom_errorbar(
    aes(
      ymin = pmax(0, mean_detected - sd_detected),
      ymax = mean_detected + sd_detected
    ),
    position = position_dodge(width = 0.78),
    width = 0.18,
    linewidth = 0.45
  ) +
  facet_wrap(
    vars(species),
    nrow = 1,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = method_colors
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(
      mult = c(0, 0.08)
    )
  ) +
  labs(
    x = NULL,
    y = "Detected genes per library",
    fill = "Quantification method",
    title = "Gene detection by quantification method",
    subtitle = "Bars show mean ± SD across libraries"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(
      fill = "#E5E5E5",
      colour = "black"
    ),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    )
  )

# -------------------------------------------------------------------------
# Shared and method-specific genes
# A group-level feature must be detected in at least half the libraries:
# 2 of 3 libraries for the present design.
# -------------------------------------------------------------------------

group_sizes <- metadata %>%
  count(
    species,
    condition,
    name = "group_size"
  ) %>%
  mutate(
    minimum_replicates = ceiling(group_size / 2)
  )

group_detection <- detection %>%
  group_by(
    species,
    condition,
    method,
    gene_id
  ) %>%
  summarise(
    detected_replicates = sum(detected),
    .groups = "drop"
  ) %>%
  left_join(
    group_sizes,
    by = c(
      "species",
      "condition"
    )
  ) %>%
  mutate(
    group_detected = (
      detected_replicates >= minimum_replicates
    )
  ) %>%
  select(
    species,
    condition,
    method,
    gene_id,
    group_detected
  ) %>%
  pivot_wider(
    names_from = method,
    values_from = group_detected,
    values_fill = FALSE
  ) %>%
  mutate(
    detection_class = case_when(
      `STAR-featureCounts` & Salmon ~ "Shared",
      `STAR-featureCounts` & !Salmon ~
        "STAR-featureCounts only",
      !`STAR-featureCounts` & Salmon ~
        "Salmon only",
      TRUE ~ "Not detected"
    )
  ) %>%
  filter(
    detection_class != "Not detected"
  )

overlap_summary <- group_detection %>%
  count(
    species,
    condition,
    detection_class,
    name = "features"
  ) %>%
  mutate(
    detection_class = factor(
      detection_class,
      levels = c(
        "Shared",
        "STAR-featureCounts only",
        "Salmon only"
      )
    )
  )

write_tsv(
  overlap_summary,
  overlap_table_output
)

overlap_colors <- c(
  "Shared" = "#4D9221",
  "STAR-featureCounts only" = "#4B2588",
  "Salmon only" = "#D98C3F"
)

overlap_plot <- ggplot(
  overlap_summary,
  aes(
    x = condition,
    y = features,
    fill = detection_class
  )
) +
  geom_col(
    width = 0.68,
    colour = "black",
    linewidth = 0.25
  ) +
  facet_wrap(
    vars(species),
    nrow = 1,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = overlap_colors,
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(
      mult = c(0, 0.08)
    )
  ) +
  labs(
    x = NULL,
    y = "Detected genes",
    fill = "Detection class",
    title = "Agreement between quantification strategies",
    subtitle = paste(
      "Group detection requires presence in at least",
      "half of the libraries"
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(
      fill = "#E5E5E5",
      colour = "black"
    ),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    )
  )

combined_plot <- (
  method_plot /
    overlap_plot
) +
  plot_annotation(
    tag_levels = "A"
  )

ggsave(
  method_plot_output,
  method_plot,
  width = 8.5,
  height = 4.5,
  units = "in"
)

ggsave(
  overlap_plot_output,
  overlap_plot,
  width = 8.5,
  height = 4.5,
  units = "in"
)

ggsave(
  combined_plot_output,
  combined_plot,
  width = 9,
  height = 9,
  units = "in"
)