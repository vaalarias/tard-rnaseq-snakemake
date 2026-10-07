#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(VennDiagram)
    library(dplyr)
    library(grid)
    library(yaml)
})

# ============================================================
# Read inputs
# ============================================================

res <- readRDS(
    snakemake@input[["deseq_results"]]
)

config <- yaml::read_yaml(
    snakemake@input[["config"]]
)

# ============================================================
# Read plot parameters from config
# ============================================================

plot_width <- config$plots$general$width
plot_height <- config$plots$general$height
font_size <- config$plots$general$font_size

venn_colors <- unlist(
    config$plots$colors$species
)

if (length(venn_colors) < 2) {
    stop(
        "At least two colours must be defined in ",
        "config$plots$colors$species."
    )
}

species_colors <- unlist(
    config$plots$colors$species
)

required_color_names <- c(
    "gadabouti",
    "experimentalis"
)

missing_colors <- setdiff(
    required_color_names,
    names(species_colors)
)

if (length(missing_colors)) {
    stop(
        "Missing species colours in config: ",
        paste(missing_colors, collapse = ", ")
    )
}

venn_colors <- c(
    species_colors[["gadabouti"]],
    species_colors[["experimentalis"]]
)
# ============================================================
# Validate DESeq2 results
# ============================================================

required_species <- c(
    "gadabouti",
    "experimentalis"
)

missing_species <- setdiff(
    required_species,
    names(res)
)

if (length(missing_species)) {
    stop(
        "Species missing from DESeq2 results: ",
        paste(missing_species, collapse = ", ")
    )
}

required_columns <- c(
    "GeneID",
    "sig"
)

for (sp in required_species) {

    missing_columns <- setdiff(
        required_columns,
        colnames(res[[sp]])
    )

    if (length(missing_columns)) {
        stop(
            "Missing columns for ",
            sp,
            ": ",
            paste(missing_columns, collapse = ", ")
        )
    }
}

# ============================================================
# Select DEGs exactly as in the original analysis
# ============================================================

gad_res <- res[["gadabouti"]]
exp_res <- res[["experimentalis"]]

sig_gad <- gad_res %>%
    filter(
        sig == "yes"
    ) %>%
    pull(GeneID)

sig_exp <- exp_res %>%
    filter(
        sig == "yes"
    ) %>%
    pull(GeneID)

sig_list <- list(
    Gadabouti = sig_gad,
    Experimentalis = sig_exp
)

message(
    "Gadabouti DEGs supplied to VennDiagram: ",
    length(sig_gad)
)

message(
    "Experimentalis DEGs supplied to VennDiagram: ",
    length(sig_exp)
)

# ============================================================
# Create Venn diagram
# ============================================================

venn_plot <- venn.diagram(
    x = sig_list,
    filename = NULL,

    fill = venn_colors,
    alpha = 0.5,

    cex = 2,
    fontface = "bold",

    cat.cex = 1.5,
    cat.fontface = "bold",

    cat.pos = c(
        -20,
        20
    ),

    margin = 0.1
)

# ============================================================
# Save PDF
# ============================================================

output_path <- as.character(
    snakemake@output[["plot"]]
)

dir.create(
    dirname(output_path),
    recursive = TRUE,
    showWarnings = FALSE
)

pdf(
    output_path,
    width = plot_width,
    height = plot_height,
    useDingbats = FALSE
)

grid.newpage()
grid.draw(venn_plot)

dev.off()