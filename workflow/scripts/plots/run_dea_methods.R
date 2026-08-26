#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(dplyr)
  library(tibble)
  library(yaml)
})

raw <- readRDS(snakemake@input[["raw_counts"]])
cfg <- yaml::read_yaml(snakemake@input[["config"]])
cts <- raw$counts
coldata <- raw$coldata
rownames(coldata) <- coldata$sample
coldata <- coldata[colnames(cts), , drop = FALSE]
coldata$condition <- factor(coldata$condition, levels = c("control", "anhydrobiosis"))
padj_cut <- cfg$deseq2$padj_threshold
lfc_cut <- cfg$deseq2$log2fc_threshold
species <- unique(coldata$species)

prepare_dge <- function(sp) {
  idx <- coldata$species == sp
  dge <- DGEList(counts = cts[, idx, drop = FALSE], samples = coldata[idx, , drop = FALSE])
  keep <- filterByExpr(dge, group = dge$samples$condition)
  calcNormFactors(dge[keep, , keep.lib.sizes = FALSE], method = "TMM")
}

results_deseq <- results_edger <- results_limma <- list()
diagnostics <- list()

for (sp in species) {
  idx <- coldata$species == sp
  meta <- coldata[idx, , drop = FALSE]

  dds <- DESeqDataSetFromMatrix(cts[, idx, drop = FALSE], meta, design = ~ condition)
  dds <- DESeq(dds)
  ds <- as.data.frame(results(dds, contrast = c("condition", "anhydrobiosis", "control"))) %>%
    rownames_to_column("GeneID") %>%
    mutate(sig = ifelse(!is.na(padj) & padj < padj_cut &
                          abs(log2FoldChange) >= lfc_cut, "yes", "no"))
  results_deseq[[sp]] <- ds

  dge <- prepare_dge(sp)
  design <- model.matrix(~ condition, data = dge$samples)
  dge <- estimateDisp(dge, design)
  lrt <- glmLRT(glmFit(dge, design), contrast = c(0, 1))
  eg <- topTags(lrt, n = Inf, sort.by = "none")$table %>%
    rownames_to_column("GeneID") %>%
    mutate(sig = ifelse(FDR < padj_cut & abs(logFC) >= lfc_cut, "yes", "no"))
  results_edger[[sp]] <- eg

  v <- voom(dge, design, plot = FALSE)
  cm <- makeContrasts(Contraste = conditionanhydrobiosis, levels = design)
  lm <- topTable(eBayes(contrasts.fit(lmFit(v, design), cm)), coef = "Contraste",
                 n = Inf, sort.by = "none") %>%
    rownames_to_column("GeneID") %>%
    mutate(sig = ifelse(adj.P.Val < padj_cut & abs(logFC) >= lfc_cut, "yes", "no"))
  results_limma[[sp]] <- lm
  diagnostics[[sp]] <- list(dge = dge, design = design)
}

saveRDS(results_deseq, snakemake@output[["deseq2"]])
saveRDS(results_edger, snakemake@output[["edger"]])
saveRDS(results_limma, snakemake@output[["limma"]])
saveRDS(diagnostics, snakemake@output[["diagnostics"]])

get_sig <- function(x, method) {
  unique(x$GeneID[!is.na(x$sig) & x$sig == "yes"])
}
summary <- bind_rows(lapply(species, function(sp) {
  sets <- list(DESeq2 = get_sig(results_deseq[[sp]]),
               edgeR = get_sig(results_edger[[sp]]),
               limma_voom = get_sig(results_limma[[sp]]))
  data.frame(
    species = sp,
    DESeq2 = length(sets$DESeq2), edgeR = length(sets$edgeR),
    limma_voom = length(sets$limma_voom),
    shared_all = length(Reduce(intersect, sets)),
    DESeq2_edgeR = length(intersect(sets$DESeq2, sets$edgeR)),
    DESeq2_limma = length(intersect(sets$DESeq2, sets$limma_voom)),
    edgeR_limma = length(intersect(sets$edgeR, sets$limma_voom))
  )
}))
write.csv(summary, snakemake@output[["summary"]], row.names = FALSE)

