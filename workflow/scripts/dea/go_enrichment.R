#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ggplot2)
  library(yaml)
})


# ------------------------------
#   Inputs Snakemake
# ------------------------------
annot_exp  <- readRDS(snakemake@input$annot_exp)
annot_gad  <- readRDS(snakemake@input$annot_gad)
term2gene  <- readRDS(snakemake@input$term2gene)
term2name  <- readRDS(snakemake@input$term2name)
config     <- yaml::read_yaml(snakemake@input$config)

pdf_exp_up   <- snakemake@output$exp_up_pdf
pdf_exp_down <- snakemake@output$exp_down_pdf
pdf_gad_up   <- snakemake@output$gad_up_pdf
pdf_gad_down <- snakemake@output$gad_down_pdf

# ------------------------------
#   Funciones
# ------------------------------

clean_gene_id <- function(x){
  x <- gsub("^LOC", "", x)      # quitar LOC
  x <- trimws(x)                # quitar espacios
  toupper(x)
}

annot_exp$GeneID <- clean_gene_id(annot_exp$GeneID)
annot_gad$GeneID <- clean_gene_id(annot_gad$GeneID)
term2gene$Gene <- clean_gene_id(term2gene$Gene)

# Separar up/down
padj_thr <- config$plots$thresholds$padj
lfc_thr  <- config$plots$thresholds$log2fc

annot_exp$regulation <- "nonsig"
annot_exp$regulation[annot_exp$padj < padj_thr & annot_exp$log2FoldChange > lfc_thr] <- "up"
annot_exp$regulation[annot_exp$padj < padj_thr & annot_exp$log2FoldChange < -lfc_thr] <- "down"

annot_gad$regulation <- "nonsig"
annot_gad$regulation[annot_gad$padj < padj_thr & annot_gad$log2FoldChange > lfc_thr] <- "up"
annot_gad$regulation[annot_gad$padj < padj_thr & annot_gad$log2FoldChange < -lfc_thr] <- "down"

# Función segura para dotplot
run_enrich <- function(genes, species_name, outfile) {

  dir.create(
    dirname(outfile),
    recursive = TRUE,
    showWarnings = FALSE
  )

  genes <- unique(na.omit(genes))

  ego <- tryCatch(
    enricher(
      gene = genes,
      TERM2GENE = term2gene,
      TERM2NAME = term2name
    ),
    error = function(e) {
      message(
        "Enricher failed for ",
        species_name,
        ": ",
        e$message
      )
      NULL
    }
  )

  no_results <- (
    is.null(ego) ||
    nrow(ego@result) == 0 ||
    all(is.na(ego@result$ID)) ||
    all(ego@result$ID == "")
  )

  pdf(
    outfile,
    width = config$plots$general$width,
    height = config$plots$general$height
  )

  if (no_results) {

    message("No enrichment results for ", species_name)

    plot.new()
    title(
      main = paste(
        "No enrichment results for",
        species_name
      )
    )

  } else {

    print(
      dotplot(
        ego,
        showCategory = config$plots$general$showCategory,
        title = paste("GO enrichment:", species_name)
      )
    )
  }

  dev.off()
}

# ------------------------------
#   Ejecutar para cada conjunto
# ------------------------------

run_enrich(
  genes = annot_exp$GeneID[annot_exp$regulation=="up"],
  species_name = "experimentalis_up",
  outfile = pdf_exp_up
)

run_enrich(
  genes = annot_exp$GeneID[annot_exp$regulation=="down"],
  species_name = "experimentalis_down",
  outfile = pdf_exp_down
)

run_enrich(
  genes = annot_gad$GeneID[annot_gad$regulation=="up"],
  species_name = "gadabouti_up",
  outfile = pdf_gad_up
)

run_enrich(
  genes = annot_gad$GeneID[annot_gad$regulation=="down"],
  species_name = "gadabouti_down",
  outfile = pdf_gad_down
)
