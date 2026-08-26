#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

colors_condition <- c(control = "#1b9e77", anhydrobiosis = "#d95f02", rehydrated = "#d95f02")
histone_cols <- c(H1 = "#318bc8ff", H2A = "#277d21ff", H2B = "#ad181bff",
                  H3 = "#f78009ff", H4 = "#6A3D9A", Other = "gray70")
epigenetic_cols <- c(Polycomb = "#6A3D9A", KMT = "#277d21ff", KDM = "#D95F02",
  HAT = "#7570B3", HDAC = "#E7298A", Epitranscriptomic = "#009E73", Other = "gray70")
deg_cols <- c(DEG = "#d73027", not_DEG = "gray90")

raw <- readRDS(snakemake@input[["raw_counts"]])
coldata <- raw$coldata
rownames(coldata) <- coldata$sample
tpm <- as.matrix(readRDS(snakemake@input[["tpm"]]))
vst_obj <- readRDS(snakemake@input[["vst"]])
vst <- if (is.matrix(vst_obj)) {
  vst_obj
} else {
  as.matrix(
    SummarizedExperiment::assay(vst_obj)
  )
}
res <- readRDS(snakemake@input[["deseq_results"]])
out_hist <- snakemake@output[["histone_dir"]]
out_epi <- snakemake@output[["epigenetic_dir"]]
dir.create(out_hist, recursive = TRUE, showWarnings = FALSE)
dir.create(out_epi, recursive = TRUE, showWarnings = FALSE)

clean_label <- function(x) {
  x <- gsub("(?i)histone[_ ]*", "", x, perl = TRUE)
  x <- gsub("(?i) protein", "", x, perl = TRUE)
  x <- gsub("(?i)putative ", "", x, perl = TRUE)
  trimws(gsub("_", " ", x))
}
hist_family <- function(x) case_when(
  grepl("H2A", x, ignore.case = TRUE) ~ "H2A", grepl("H2B", x, ignore.case = TRUE) ~ "H2B",
  grepl("H1", x, ignore.case = TRUE) ~ "H1", grepl("H3", x, ignore.case = TRUE) ~ "H3",
  grepl("H4", x, ignore.case = TRUE) ~ "H4", TRUE ~ "Other")
epi_family <- function(x) { p <- tolower(x); case_when(
  grepl("polycomb|prc1|prc2|ezh|suz12|eed|ring1|bmi1", p) ~ "Polycomb",
  grepl("demethylase|jumonji|kdm", p) ~ "KDM",
  grepl("methyltransferase|set domain|kmt|smyd", p) ~ "KMT",
  grepl("deacetylase|hdac|sirtuin", p) ~ "HDAC",
  grepl("acetyltransferase|gnat|myst|hat", p) ~ "HAT", TRUE ~ "Other") }
load_features <- function(path, family_fun) {
  x <- read.delim(path, header = FALSE, stringsAsFactors = FALSE)[, 1:2, drop = FALSE]
  colnames(x) <- c("gene_id", "tag")
  x <- x[!duplicated(x$gene_id) & x$gene_id %in% rownames(tpm), , drop = FALSE]
  x$Family <- family_fun(x$tag); x$label <- clean_label(x$tag); x
}
hist <- load_features(snakemake@input[["histones"]], hist_family)
histone_order <- c(
  "H1",
  "H2A",
  "H2B",
  "H3",
  "H4",
  "Other"
)

hist$Family <- factor(
  hist$Family,
  levels = histone_order
)

hist <- hist %>%
  arrange(
    Family,
    label,
    gene_id
  )

epi <- load_features(snakemake@input[["epigenetic_regulators"]], epi_family)
epi_rna_raw <- read.delim(
  snakemake@input[["epitranscriptomic"]],
  header = FALSE,
  stringsAsFactors = FALSE,
  fill = TRUE,
  quote = ""
)

if (ncol(epi_rna_raw) < 1) {
  stop("The epitranscriptomic table is empty.")
}

# Soporta temporalmente archivos viejos con solamente una columna
if (ncol(epi_rna_raw) == 1) {
  colnames(epi_rna_raw) <- "gene_id"

  epi_rna_raw$tag <- NA_character_
} else {
  epi_rna_raw <- epi_rna_raw[, 1:2, drop = FALSE]

  colnames(epi_rna_raw) <- c(
    "gene_id",
    "tag"
  )
}

epi_rna_raw$gene_id <- trimws(
  as.character(epi_rna_raw$gene_id)
)

epi_rna_raw$tag <- trimws(
  as.character(epi_rna_raw$tag)
)

# Recover label
existing_labels <- epi %>%
  select(
    gene_id,
    existing_tag = tag
  ) %>%
  distinct(
    gene_id,
    .keep_all = TRUE
  )

epi_rna <- epi_rna_raw %>%
  left_join(
    existing_labels,
    by = "gene_id"
  ) %>%
  mutate(
    tag = case_when(
      !is.na(tag) & nzchar(tag) ~ tag,
      !is.na(existing_tag) & nzchar(existing_tag) ~ existing_tag,
      TRUE ~ gene_id
    ),

    Family = "Epitranscriptomic",

    label = case_when(
      tag == gene_id ~ gene_id,
      TRUE ~ clean_label(tag)
    )
  ) %>%
  filter(
    gene_id %in% rownames(tpm)
  ) %>%
  select(
    gene_id,
    tag,
    Family,
    label
  ) %>%
  distinct(
    gene_id,
    .keep_all = TRUE
  )

epi_ids <- epi_rna$gene_id

epi <- bind_rows(
  epi[!epi$gene_id %in% epi_ids, , drop = FALSE],
  epi_rna
)

deg_ids <- lapply(res, function(x) unique(x$GeneID[!is.na(x$sig) & x$sig == "yes"]))
row_z <- function(x) { z <- t(scale(t(x))); z[is.na(z)] <- 0; z }
col_tpm <- function(x) {
  q <- quantile(as.numeric(x), c(0, .5, .8, .95, 1), na.rm = TRUE)
  if (length(unique(q)) < 5) q <- seq(min(x), max(x) + 1e-9, length.out = 5)
  colorRamp2(q, c("white", "#ffffbf", "#fdae61", "#d73027", "#67001f"))
}
col_z <- function(x) { lim <- min(max(abs(x), na.rm = TRUE), 2); if (!is.finite(lim) || lim == 0) lim <- 2
  colorRamp2(c(-lim, 0, lim), c("#313695", "white", "#a50026")) }

collapse_median <- function(mat, ann, summary_file) {
  groups <- unique(ifelse(is.na(ann$label) | !nzchar(ann$label), ann$gene_id, ann$label))
  labels <- ifelse(is.na(ann$label) | !nzchar(ann$label), ann$gene_id, ann$label)
  out <- do.call(rbind, lapply(groups, function(g) apply(mat[labels == g, , drop = FALSE], 2, median)))
  rownames(out) <- groups
  ann2 <- bind_rows(lapply(groups, function(g) {
    i <- which(labels == g)
    data.frame(gene_id = g, tag = g, Family = ann$Family[i[1]], label = g,
      DEG = ifelse(any(ann$DEG[i] == "DEG"), "DEG", "not_DEG"),
      LOCs = paste(ann$gene_id[i], collapse = ";"),
      DEG_LOCs = paste(ann$gene_id[i][ann$DEG[i] == "DEG"], collapse = ";"))
  }))
  write.csv(transmute(ann2, Label = label, `LOCs in this label` = LOCs,
                      `LOCs differentially expressed in this label` = DEG_LOCs),
            summary_file, row.names = FALSE)
  list(mat = out, ann = ann2)
}

plot_feature <- function(expr, ann, sp, mode, colors, outfile, title,
                         collapse = FALSE, summary_file = NULL,
                         cluster_rows = FALSE, width = 8, height = 10) {
  samples <- intersect(coldata$sample[coldata$species == sp], colnames(expr))
  genes <- intersect(ann$gene_id, rownames(expr))
  a <- ann[match(genes, ann$gene_id), , drop = FALSE]
  m <- expr[a$gene_id, samples, drop = FALSE]
  if (mode != "zscore_VST") { keep <- rowSums(m > 0) > 0; m <- m[keep,,drop=FALSE]; a <- a[keep,,drop=FALSE] }
  a$DEG <- ifelse(a$gene_id %in% deg_ids[[sp]], "DEG", "not_DEG")
  if (collapse) { x <- collapse_median(m, a, summary_file); m <- x$mat; a <- x$ann }
  a$Family <- factor(
  a$Family,
  levels = names(colors)
  )

  if (!cluster_rows) {
    row_order <- order(
      a$Family,
      a$label,
      a$gene_id
    )

    m <- m[
      row_order,
      ,
      drop = FALSE
    ]

    a <- a[
      row_order,
      ,
      drop = FALSE
    ]
  }
  if (!nrow(m)) return(invisible(NULL))
  if (mode == "logTPM") { pm <- log10(m + 1); cf <- col_tpm(pm); lname <- "logTPM" }
  else if (mode == "zscore_logTPM") { pm <- row_z(log10(m + 1)); cf <- col_z(pm); lname <- "zTPM" }
  else { pm <- row_z(m); cf <- col_z(pm); lname <- "zVST" }
  top <- HeatmapAnnotation(Condition = coldata$condition[match(colnames(pm), coldata$sample)],
                           col = list(Condition = colors_condition))
  left <- rowAnnotation(Family = a$Family, DEG = a$DEG,
                        col = list(Family = colors, DEG = deg_cols), na_col = "gray90")
  ht <- Heatmap(pm, name = lname, col = cf, top_annotation = top, left_annotation = left,
    row_labels = a$label, show_row_names = TRUE, show_column_names = TRUE,
    cluster_rows = cluster_rows, cluster_columns = TRUE,
    column_title = paste(title, sp, mode, sep = " | "),
    row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 10))
  pdf(outfile, width = width, height = height); grid.newpage(); draw(ht); dev.off()
}

species <- c("gadabouti", "experimentalis")
for (sp in species) {
  for (mode in c("logTPM", "zscore_logTPM", "zscore_VST")) {
    expr <- if (mode == "zscore_VST") vst else tpm
    plot_feature(expr, hist, sp, mode, histone_cols,
      file.path(out_hist, paste0("histones_", sp, "_", mode, "_LOCs.pdf")), "Histones LOC-level")
    plot_feature(expr, hist, sp, mode, histone_cols,
      file.path(out_hist, paste0("histones_", sp, "_", mode, "_medianCollapsed.pdf")),
      "Histones median-collapsed", TRUE,
      file.path(out_hist, paste0("histones_", sp, "_collapsed_summary.csv")), width = 7, height = 8)
  }

  for (cat in c("Polycomb", "KMT", "KDM", "HAT", "HDAC", "Epitranscriptomic")) {
    a <- epi[epi$Family == cat, , drop = FALSE]
    if (!nrow(a)) next
    for (mode in c("logTPM", "zscore_logTPM", "zscore_VST")) {
      expr <- if (mode == "zscore_VST") vst else tpm
      plot_feature(expr, a, sp, mode, epigenetic_cols,
        file.path(out_epi, paste0("epigenetic_", cat, "_", sp, "_", mode, "_LOCs.pdf")),
        paste("Epigenetic", cat, "LOC-level"), width = 8, height = 8)
      if (mode != "zscore_VST") plot_feature(expr, a, sp, mode, epigenetic_cols,
        file.path(out_epi, paste0("epigenetic_", cat, "_", sp, "_", mode, "_medianCollapsed.pdf")),
        paste("Epigenetic", cat, "median-collapsed"), TRUE,
        file.path(out_epi, paste0("epigenetic_", cat, "_", sp, "_collapsed_summary.csv")),
        width = 8, height = 8)
    }
  }

  for (mode in c("logTPM", "zscore_logTPM", "zscore_VST")) {
    expr <- if (mode == "zscore_VST") vst else tpm
    plot_feature(expr, epi, sp, mode, epigenetic_cols,
      file.path(out_epi, paste0("epigenetic_ALL_", sp, "_", mode, "_LOCs.pdf")),
      "Epigenetic all categories LOC-level", width = 9, height = 14)
    plot_feature(expr, epi, sp, mode, epigenetic_cols,
      file.path(out_epi, paste0("epigenetic_ALL_", sp, "_", mode, "_medianCollapsed.pdf")),
      "Epigenetic all categories median-collapsed", TRUE,
      file.path(out_epi, paste0("epigenetic_ALL_", sp, "_collapsed_summary.csv")), width = 9, height = 12)
    if (mode != "logTPM") plot_feature(expr, epi, sp, mode, epigenetic_cols,
      file.path(out_epi, paste0("epigenetic_ALL_", sp, "_", mode, "_LOCs_rowClustered.pdf")),
      "Epigenetic all categories LOC-level row-clustered", cluster_rows = TRUE,
      width = 9, height = 14)
  }
}

