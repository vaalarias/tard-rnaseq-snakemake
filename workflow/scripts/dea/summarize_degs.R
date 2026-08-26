#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(yaml)
})

### Inputs

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

padj_threshold <- config$deseq2$padj_threshold
log2fc_threshold <- config$deseq2$log2fc_threshold


### Validate columns

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
      paste(missing_columns, collapse = ", ")
    )
  }
}

validate_results(exp_results, "experimentalis")
validate_results(gad_results, "gadabouti")


### Classify differential expression

classify_degs <- function(x) {

  x %>%
    mutate(
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

exp_results <- classify_degs(exp_results)
gad_results <- classify_degs(gad_results)


### Generate gene sets

gene_sets <- list(
  experimentalis = list(
    up = unique(
      exp_results$GeneID[
        exp_results$regulation == "up"
      ]
    ),
    down = unique(
      exp_results$GeneID[
        exp_results$regulation == "down"
      ]
    )
  ),

  gadabouti = list(
    up = unique(
      gad_results$GeneID[
        gad_results$regulation == "up"
      ]
    ),
    down = unique(
      gad_results$GeneID[
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


### Shared genes

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


### Summary by species

species_summary <- data.frame(
  species = c(
    "experimentalis",
    "gadabouti"
  ),
  upregulated = c(
    length(gene_sets$experimentalis$up),
    length(gene_sets$gadabouti$up)
  ),
  downregulated = c(
    length(gene_sets$experimentalis$down),
    length(gene_sets$gadabouti$down)
  )
)

species_summary$total_degs <- (
  species_summary$upregulated +
  species_summary$downregulated
)


### Shared summary

shared_summary <- data.frame(
  category = c(
    "Shared upregulated",
    "Shared downregulated",
    "Shared same direction",
    "Shared opposite direction",
    "Experimentalis up / Gadabouti down",
    "Experimentalis down / Gadabouti up",
    "Shared DEGs regardless of direction"
  ),
  number_of_genes = c(
    length(gene_sets$shared$up),
    length(gene_sets$shared$down),
    length(gene_sets$shared$same_direction),
    length(gene_sets$shared$opposite_direction),
    length(gene_sets$shared$exp_up_gad_down),
    length(gene_sets$shared$exp_down_gad_up),
    length(gene_sets$shared$all)
  )
)


### Gene-level comparison table

all_deg_ids <- union(
  gene_sets$experimentalis$all,
  gene_sets$gadabouti$all
)

gene_comparison <- data.frame(
  GeneID = all_deg_ids
) %>%
  left_join(
    exp_results %>%
      select(
        GeneID,
        experimentalis_log2FC = log2FoldChange,
        experimentalis_padj = padj,
        experimentalis_regulation = regulation
      ),
    by = "GeneID"
  ) %>%
  left_join(
    gad_results %>%
      select(
        GeneID,
        gadabouti_log2FC = log2FoldChange,
        gadabouti_padj = padj,
        gadabouti_regulation = regulation
      ),
    by = "GeneID"
  ) %>%
  arrange(GeneID)


### Save outputs

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

message("\nDifferential expression summary:")
print(species_summary)

message("\nShared DEG summary:")
print(shared_summary)