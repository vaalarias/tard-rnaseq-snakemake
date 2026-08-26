#!/usr/bin/env Rscript

library(DESeq2)
library(tidyverse)
library(yaml)

### ---------------------------
### 1. Load config and raw counts
### ---------------------------

config <- yaml::read_yaml(snakemake@input$config)
contrast_list   <- config$deseq2$contrasts[[1]]
padj_threshold  <- config$deseq2$padj_threshold
log2fc_threshold <- config$deseq2$log2fc_threshold

raw <- readRDS(snakemake@input$raw_counts)
cts <- raw$counts
coldata <- raw$coldata

### ---------------------------
### 2. Build DESeq2 object
### ---------------------------

dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design = as.formula(config$deseq2$design)
)

### ---------------------------
### 3. Split by species
### ---------------------------

species_list <- unique(coldata$species)
dds_list <- lapply(species_list, function(sp){
  dds_sp <- dds[, dds$species == sp]
  dds_sp$species <- droplevels(dds_sp$species)
  dds_sp$condition <- droplevels(dds_sp$condition)
  dds_sp
})
names(dds_list) <- species_list

### ---------------------------
### 4. Run DESeq2 per species
### ---------------------------

results_all <- list()

for(sp in names(dds_list)){
  message("Running DESeq2 for species: ", sp)
  dds_sp <- dds_list[[sp]]
  
  # Overwrite design for condition only
  design(dds_sp) <- ~ condition
  dds_sp <- DESeq(dds_sp)
  
  res <- results(dds_sp, contrast = unlist(contrast_list))
  res_df <- as.data.frame(res) %>% rownames_to_column("GeneID")
  
  res_df <- res_df %>%
    mutate(sig = ifelse(padj < padj_threshold & abs(log2FoldChange) >= log2fc_threshold, "yes", "no"))
  
  results_all[[sp]] <- res_df
}

### ---------------------------
### 5. Save results
### ---------------------------

saveRDS(results_all, snakemake@output$deseq_results)

out_dir <- snakemake@params$outdir
for(sp in names(results_all)){
  write.csv(results_all[[sp]], file.path(out_dir, paste0(sp,"_deseq2_results.csv")), row.names=FALSE)
}
