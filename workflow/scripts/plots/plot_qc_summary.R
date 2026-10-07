#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
})

qc <- read.csv(snakemake@input[["qc_summary"]])

colors_condition <- c(
    control = "#1b9e77",
    anhydrobiosis = "#d95f02"
)

set.seed(123)

p1 <- ggplot(
    qc,
    aes(x = species, y = million_reads, fill = condition)
) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
        aes(color = condition),
        width = 0.15,
        size = 2.5
    ) +
    scale_fill_manual(values = colors_condition) +
    scale_color_manual(values = colors_condition) +
    labs(
        x = NULL,
        y = "Reads (millions)",
        fill = "Condition",
        color = "Condition"
    ) +
    theme_minimal(base_size = 12)

p2 <- ggplot(
    qc,
    aes(x = species, y = features_above_threshold, fill = condition)
) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
        aes(color = condition),
        width = 0.15,
        size = 2.5
    ) +
    scale_fill_manual(values = colors_condition) +
    scale_color_manual(values = colors_condition) +
    labs(
        x = NULL,
        y = "Detected Features > 10 Reads",
        fill = "Condition",
        color = "Condition"
    ) +
    theme_minimal(base_size = 12)

p3 <- ggplot(
    qc,
    aes(x = species, y = percent_aligned, fill = condition)
) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
        aes(color = condition),
        width = 0.15,
        size = 2.5
    ) +
    scale_fill_manual(values = colors_condition) +
    scale_color_manual(values = colors_condition) +
    labs(
        x = NULL,
        y = "Overall alignment rate (%)",
        fill = "Condition",
        color = "Condition"
    ) +
    theme_minimal(base_size = 12)

qc_plot <- (
    p3 + p1 + p2 +
        plot_layout(ncol = 3, guides = "collect") +
        plot_annotation(tag_levels = "A")
) &
    theme(
        legend.position = "bottom",
        plot.tag = element_text(face = "bold", size = 14)
    )

ggsave(
    filename = snakemake@output[["plot"]],
    plot = qc_plot,
    width = 18,
    height = 5,
    dpi = 300
)