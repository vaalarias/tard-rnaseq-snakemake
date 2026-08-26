#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

vst_list <- readRDS(snakemake@input[["vst_list"]])
tpm <- as.matrix(readRDS(snakemake@input[["tpm"]]))
raw <- readRDS(snakemake@input[["raw_counts"]])
coldata <- raw$coldata
colors <- c(control = "#1b9e77", anhydrobiosis = "#d95f02")
vst <- lapply(vst_list, assay)

plot_comparison <- function(mat_exp, mat_gad, meta_exp, meta_gad, condition,
                            title, outfile) {
  genes <- intersect(rownames(mat_exp), rownames(mat_gad))
  df <- data.frame(
    experimentalis = rowMeans(mat_exp[genes, meta_exp$condition == condition, drop = FALSE]),
    gadabouti = rowMeans(mat_gad[genes, meta_gad$condition == condition, drop = FALSE])
  )
  r <- cor(df$experimentalis, df$gadabouti, method = "pearson")
  p <- ggplot(df, aes(experimentalis, gadabouti)) +
    geom_point(color = colors[[condition]], alpha = .30, size = .4) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    coord_equal() +
    annotate("text", x = Inf, y = -Inf, label = paste0("Pearson r = ", round(r, 3)),
             hjust = 1.1, vjust = -.5, size = 3) +
    labs(title = title, x = "experimentalis", y = "gadabouti") +
    theme_classic(base_size = 12)
  ggsave(outfile, p, width = 6, height = 6)
}

meta_exp <- as.data.frame(colData(vst_list[["experimentalis"]]))
meta_gad <- as.data.frame(colData(vst_list[["gadabouti"]]))
tpm_exp <- log2(tpm[, coldata$species == "experimentalis", drop = FALSE] + 1)
tpm_gad <- log2(tpm[, coldata$species == "gadabouti", drop = FALSE] + 1)
meta_exp_tpm <- coldata[coldata$species == "experimentalis", , drop = FALSE]
meta_gad_tpm <- coldata[coldata$species == "gadabouti", , drop = FALSE]

plot_comparison(vst$experimentalis, vst$gadabouti, meta_exp, meta_gad, "control",
                "Control vs Control (VST)", snakemake@output[["control_vst"]])
plot_comparison(tpm_exp, tpm_gad, meta_exp_tpm, meta_gad_tpm, "control",
                "Control vs Control (log2 TPM + 1)", snakemake@output[["control_tpm"]])
plot_comparison(vst$experimentalis, vst$gadabouti, meta_exp, meta_gad, "anhydrobiosis",
                "Anhydrobiosis vs Anhydrobiosis (VST)", snakemake@output[["anh_vst"]])
plot_comparison(tpm_exp, tpm_gad, meta_exp_tpm, meta_gad_tpm, "anhydrobiosis",
                "Anhydrobiosis vs Anhydrobiosis (log2 TPM + 1)", snakemake@output[["anh_tpm"]])

