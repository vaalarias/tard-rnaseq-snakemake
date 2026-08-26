#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tibble)
})

colors_condition <- c(control = "#1b9e77", anhydrobiosis = "#d95f02")
colors_species <- c(experimentalis = "#56b4e9", gadabouti = "#4b2588ff")
histone_cols <- c(H1 = "#318bc8ff", H2A = "#277d21ff", H2B = "#ad181bff",
                  H3 = "#f78009ff", H4 = "#6A3D9A", Other = "gray70")
family_order <- names(histone_cols)

tpm <- as.matrix(readRDS(snakemake@input[["tpm"]]))
raw <- readRDS(snakemake@input[["raw_counts"]])
coldata <- raw$coldata
hist <- read.delim(snakemake@input[["histones"]], header = FALSE,
                   stringsAsFactors = FALSE)[, 1:2, drop = FALSE]
colnames(hist) <- c("gene_id", "tag")
hist <- hist[!duplicated(hist$gene_id) & hist$gene_id %in% rownames(tpm), , drop = FALSE]
if (!nrow(hist)) stop("No histone genes overlap the TPM matrix")

clean_label <- function(x) gsub("(?i)histone[_ ]*", "", x, perl = TRUE)
family_of <- function(x) case_when(
  grepl("H1", x, ignore.case = TRUE) ~ "H1",
  grepl("H2A", x, ignore.case = TRUE) ~ "H2A",
  grepl("H2B", x, ignore.case = TRUE) ~ "H2B",
  grepl("H3", x, ignore.case = TRUE) ~ "H3",
  grepl("H4", x, ignore.case = TRUE) ~ "H4",
  TRUE ~ "Other"
)
zscore_rows <- function(x) {
  z <- t(scale(t(x)))
  z[is.na(z)] <- 0
  z
}
make_annotations <- function(mat, labels) {
  fam <- factor(family_of(labels), levels = family_order)
  ord <- order(as.numeric(fam), labels)
  ann_col <- coldata[match(colnames(mat), coldata$sample), , drop = FALSE]
  list(
    order = ord,
    row = rowAnnotation(Family = fam[ord], col = list(Family = histone_cols),
                        annotation_name_gp = gpar(fontsize = 9), na_col = "gray90"),
    top = HeatmapAnnotation(Condition = ann_col$condition,
                            col = list(Condition = colors_condition),
                            annotation_name_gp = gpar(fontsize = 10))
  )
}
draw_histone <- function(mat, labels, outfile, title, name, col_fun,
                         cluster_columns = TRUE, show_rows = TRUE) {
  a <- make_annotations(mat, labels)
  mat <- mat[a$order, , drop = FALSE]
  labels <- labels[a$order]
  ht <- Heatmap(mat, name = name, col = col_fun, top_annotation = a$top,
                left_annotation = a$row, row_labels = labels,
                show_row_names = show_rows, show_column_names = TRUE,
                cluster_rows = FALSE, cluster_columns = cluster_columns,
                column_title = title, column_names_gp = gpar(fontsize = 10))
  pdf(outfile, width = 8, height = 8)
  draw(ht)
  dev.off()
}

species <- c("experimentalis", "gadabouti")
z_uncollapsed_out <- setNames(as.character(snakemake@output[["z_uncollapsed"]]), species)
tpm_out <- setNames(as.character(snakemake@output[["tpm_species"]]), species)
logtpm_out <- setNames(as.character(snakemake@output[["logtpm_species"]]), species)
z_merged_out <- setNames(as.character(snakemake@output[["z_merged"]]), species)
polya_no_out <- setNames(as.character(snakemake@output[["polya_no"]]), species)
polya_yes_out <- setNames(as.character(snakemake@output[["polya_yes"]]), species)

for (sp in species) {
  sample_ids <- coldata$sample[coldata$species == sp]
  genes <- intersect(hist$gene_id, rownames(tpm))
  mat <- tpm[genes, sample_ids, drop = FALSE]
  tags <- hist$tag[match(genes, hist$gene_id)]
  labels <- clean_label(tags)

  draw_histone(zscore_rows(mat), labels, z_uncollapsed_out[[sp]],
               paste0("Histone heatmap (", sp, ")"), "Z-score",
               colorRamp2(c(-4, 0, 4), c("#313695", "white", "#a50026")))

  logmat <- log10(mat + 1)
  keep <- apply(logmat, 1, sd) > 0
  q <- quantile(logmat[keep, , drop = FALSE], c(.05, .5, .95), na.rm = TRUE)
  draw_histone(logmat[keep, , drop = FALSE], labels[keep], tpm_out[[sp]],
               paste0("Histone TPM heatmap (", sp, ")"), "tpm",
               colorRamp2(q, c("white", "#fdae61", "#a50026")))

  draw_histone(log2(mat + 1), tags, logtpm_out[[sp]],
               paste0("Histones TPM (", sp, ")"), "log2(TPM+1)",
               colorRamp2(range(log2(mat + 1)), c("white", "#a50026")))

  # Final notebook version: collapse every repeated clean label by median TPM.
  collapsed <- as.data.frame(mat) %>% rownames_to_column("gene_id") %>%
    mutate(Label = labels) %>% group_by(Label) %>%
    summarise(across(where(is.numeric), median), .groups = "drop")
  merged_labels <- collapsed$Label
  merged <- as.matrix(dplyr::select(collapsed, -Label))
  rownames(merged) <- merged_labels
  draw_histone(zscore_rows(merged), merged_labels, z_merged_out[[sp]],
               paste0("Histone heatmap (", sp, ")"), "Z-score",
               colorRamp2(c(-4, 0, 4), c("#313695", "white", "#a50026")))

  polya <- grepl("H3\\.3|H2A\\.Z|H1A|H2A-beta", tags, ignore.case = TRUE)
  for (status in c(FALSE, TRUE)) {
    idx <- which(polya == status)
    if (!length(idx)) stop("No histones for PolyA status ", status, " in ", sp)
    x <- as.data.frame(mat[idx, , drop = FALSE]) %>% rownames_to_column("gene_id") %>%
      mutate(Label = labels[idx]) %>% group_by(Label) %>%
      summarise(across(where(is.numeric), median), .groups = "drop")
    labs <- x$Label
    m <- as.matrix(dplyr::select(x, -Label)); rownames(m) <- labs
    outfile <- if (status) polya_yes_out[[sp]] else polya_no_out[[sp]]
    title <- if (status) paste0("PolyA variants (", sp, ")") else paste0("Replicative histones (", sp, ")")
    draw_histone(zscore_rows(m), labs, outfile, title, "Z-score",
                 colorRamp2(c(-4, 0, 4), c("#313695", "white", "#a50026")))
  }
}

# All-sample notebook heatmaps.
genes <- hist$gene_id
mat <- tpm[genes, , drop = FALSE]
labels <- genes
fam <- factor(family_of(hist$tag), levels = family_order)
ord <- order(fam)
ann_col <- coldata[match(colnames(mat), coldata$sample), , drop = FALSE]
top <- HeatmapAnnotation(Species = ann_col$species, Condition = ann_col$condition,
  col = list(Species = colors_species, Condition = colors_condition))
left <- rowAnnotation(Family = fam[ord], col = list(Family = histone_cols))

ht <- Heatmap(mat[ord, , drop = FALSE], name = "TPM",
  col = colorRamp2(range(mat), c("white", "red")), top_annotation = top,
  left_annotation = left, row_labels = labels[ord], show_row_names = FALSE,
  show_column_names = TRUE, cluster_rows = FALSE, cluster_columns = FALSE,
  column_title = paste0("Histone Heatmap -- ", nrow(mat), " predicted genes"))
pdf(snakemake@output[["tpm_all"]], width = 8, height = 8); draw(ht); dev.off()

z <- zscore_rows(mat)
ht <- Heatmap(z[ord, , drop = FALSE], name = "Z-score",
  col = colorRamp2(c(-4, 0, 4), c("#313695", "white", "#a50026")),
  top_annotation = top, left_annotation = left, row_labels = labels[ord],
  show_row_names = TRUE, show_column_names = TRUE, cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_title = paste0("Histone Heatmap -- ", nrow(z), " predicted genes"))
pdf(snakemake@output[["z_all"]], width = 8, height = 8); draw(ht); dev.off()

