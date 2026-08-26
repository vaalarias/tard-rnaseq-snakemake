#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(yaml)
})

### Inputs

annot_exp <- readRDS(snakemake@input[["annot_exp"]])
unan_exp  <- readRDS(snakemake@input[["unan_exp"]])
annot_gad <- readRDS(snakemake@input[["annot_gad"]])
unan_gad  <- readRDS(snakemake@input[["unan_gad"]])

config <- yaml::read_yaml(
  snakemake@input[["config"]]
)

padj_threshold <- config$deseq2$padj_threshold
log2fc_threshold <- config$deseq2$log2fc_threshold


### Classify genes independently of the existing sig column

classify_degs <- function(x, annotation_status) {

  x %>%
    mutate(
      annotation_status = annotation_status,

      regulation = case_when(
        !is.na(padj) &
          padj < padj_threshold &
          log2FoldChange >= log2fc_threshold ~ "up",

        !is.na(padj) &
          padj < padj_threshold &
          log2FoldChange <= -log2fc_threshold ~ "down",

        TRUE ~ "not_significant"
      )
    )
}

exp <- bind_rows(
  classify_degs(annot_exp, "annotated"),
  classify_degs(unan_exp, "unannotated")
)

gad <- bind_rows(
  classify_degs(annot_gad, "annotated"),
  classify_degs(unan_gad, "unannotated")
)


### Keep DEGs

exp_degs <- exp %>%
  filter(regulation != "not_significant")

gad_degs <- gad %>%
  filter(regulation != "not_significant")


### Summary per species and annotation status

summarize_species <- function(x, species_name) {

  total_degs <- n_distinct(x$GeneID)

  x %>%
    group_by(annotation_status) %>%
    summarise(
      upregulated = n_distinct(GeneID[regulation == "up"]),
      downregulated = n_distinct(GeneID[regulation == "down"]),
      total_degs = n_distinct(GeneID),
      .groups = "drop"
    ) %>%
    mutate(
      species = species_name,
      percentage_of_species_degs = round(
        100 * total_degs / sum(total_degs),
        2
      )
    ) %>%
    select(
      species,
      annotation_status,
      upregulated,
      downregulated,
      total_degs,
      percentage_of_species_degs
    )
}

annotation_summary <- bind_rows(
  summarize_species(exp_degs, "experimentalis"),
  summarize_species(gad_degs, "gadabouti")
)


### Gene-level comparison

exp_status <- exp_degs %>%
  select(
    GeneID,
    exp_regulation = regulation,
    exp_annotation = annotation_status
  )

gad_status <- gad_degs %>%
  select(
    GeneID,
    gad_regulation = regulation,
    gad_annotation = annotation_status
  )

shared_details <- inner_join(
  exp_status,
  gad_status,
  by = "GeneID"
) %>%
  mutate(
    direction_comparison = case_when(
      exp_regulation == gad_regulation ~ "same_direction",
      TRUE ~ "opposite_direction"
    ),

    annotation_comparison = case_when(
      exp_annotation == "annotated" &
        gad_annotation == "annotated" ~ "annotated_in_both",

      exp_annotation == "annotated" &
        gad_annotation == "unannotated" ~
        "annotated_only_in_experimentalis",

      exp_annotation == "unannotated" &
        gad_annotation == "annotated" ~
        "annotated_only_in_gadabouti",

      TRUE ~ "unannotated_in_both"
    )
  )


### Shared summary

shared_summary <- bind_rows(
  data.frame(
    category = "Shared DEGs regardless of direction",
    number_of_genes = n_distinct(shared_details$GeneID)
  ),
  data.frame(
    category = "Shared DEGs in the same direction",
    number_of_genes = n_distinct(
      shared_details$GeneID[
        shared_details$direction_comparison == "same_direction"
      ]
    )
  ),
  data.frame(
    category = "Shared DEGs in opposite directions",
    number_of_genes = n_distinct(
      shared_details$GeneID[
        shared_details$direction_comparison == "opposite_direction"
      ]
    )
  ),
  data.frame(
    category = "Shared DEGs annotated in both",
    number_of_genes = n_distinct(
      shared_details$GeneID[
        shared_details$annotation_comparison == "annotated_in_both"
      ]
    )
  ),
  data.frame(
    category = "Shared DEGs unannotated in both",
    number_of_genes = n_distinct(
      shared_details$GeneID[
        shared_details$annotation_comparison == "unannotated_in_both"
      ]
    )
  ),
  data.frame(
    category = "Shared DEGs annotated only in experimentalis",
    number_of_genes = n_distinct(
      shared_details$GeneID[
        shared_details$annotation_comparison ==
          "annotated_only_in_experimentalis"
      ]
    )
  ),
  data.frame(
    category = "Shared DEGs annotated only in gadabouti",
    number_of_genes = n_distinct(
      shared_details$GeneID[
        shared_details$annotation_comparison ==
          "annotated_only_in_gadabouti"
      ]
    )
  )
)


### Save

write.csv(
  annotation_summary,
  snakemake@output[["annotation_summary"]],
  row.names = FALSE
)

write.csv(
  shared_summary,
  snakemake@output[["shared_summary"]],
  row.names = FALSE
)

write.csv(
  shared_details,
  snakemake@output[["shared_details"]],
  row.names = FALSE
)

message("\nAnnotation summary:")
print(annotation_summary)

message("\nShared annotation summary:")
print(shared_summary)