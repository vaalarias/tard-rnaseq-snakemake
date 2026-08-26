#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(mgsa)
    library(tibble)
})

# Snakemake inputs
deg_exp_file <- snakemake@input$exp
deg_gad_file <- snakemake@input$gad
gaf_file     <- snakemake@input$gaf

# Outputs
term2gene_file <- snakemake@output$term2gene
term2name_file <- snakemake@output$term2name
annot_exp_file <- snakemake@output$annot_exp
annot_gad_file <- snakemake@output$annot_gad
unan_exp_file  <- snakemake@output$unan_exp
unan_gad_file  <- snakemake@output$unan_gad

# Load DEGs
exp <- read.csv(deg_exp_file)
gad <- read.csv(deg_gad_file)

# Load GAF with mgsa
GAF <- readGAF(gaf_file)

# Build TERM2GENE and TERM2NAME only ONCE
sets <- GAF@sets[names(GAF@sets) != "all"]
symbols <- GAF@itemAnnotations$symbol

TERM2GENE <- do.call(
    rbind,
    lapply(names(sets), function(go){
        idx <- sets[[go]]
        if(length(idx)==0) return(NULL)
        data.frame(
            GO = go,
            Gene = symbols[idx],
            stringsAsFactors = FALSE
        )
    })
)

TERM2NAME <- data.frame(
    GO = rownames(GAF@setAnnotations)[rownames(GAF@setAnnotations)!="all"],
    Description = GAF@setAnnotations$term[rownames(GAF@setAnnotations)!="all"],
    stringsAsFactors = FALSE
)

# Save annotation tables
saveRDS(TERM2GENE, term2gene_file)
saveRDS(TERM2NAME, term2name_file)

# Annotate genes
exp$annotation <- symbols[match(exp$GeneID, symbols)]
gad$annotation <- symbols[match(gad$GeneID, symbols)]

# Separate annotated / unannotated
exp_annot <- exp %>% filter(!is.na(annotation))
exp_unan  <- exp %>% filter(is.na(annotation))

gad_annot <- gad %>% filter(!is.na(annotation))
gad_unan  <- gad %>% filter(is.na(annotation))

# Save species-specific outputs
saveRDS(exp_annot, annot_exp_file)
saveRDS(exp_unan,  unan_exp_file)
saveRDS(gad_annot, annot_gad_file)
saveRDS(gad_unan,  unan_gad_file)
