library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)

# Input / output from Snakemake

input_file <- snakemake@input[["summary"]]
output_pdf <- snakemake@output[["pdf"]]


# Read STAR trial-mapping summary


df <- read_tsv(input_file, show_col_types = FALSE)


# Add sample metadata

df <- df %>%
    mutate(
        species = case_when(
            grepl("^E", sample) ~ "experimentalis",
            grepl("^G", sample) ~ "gadabouti",
            TRUE ~ NA_character_
        ),

        condition = case_when(
            grepl("CO", sample) | grepl("CB", sample) ~ "control",
            grepl("AO", sample) | grepl("AB", sample) ~ "anhydrobiosis",
            TRUE ~ NA_character_
        )
    )

# Labels and ordering


genome_levels <- c(
    "Hyp_dujardini",
    "Hyp_exemplaris",
    "Ram_varieornatus",
    "Par_richtersi",
    "Pam_metropolitanus"
)

genome_labels <- c(
    "Hyp_dujardini"      = "H. dujardini",
    "Hyp_exemplaris"     = "H. exemplaris",
    "Ram_varieornatus"   = "R. varieornatus",
    "Par_richtersi"      = "P. richtersi",
    "Pam_metropolitanus" = "P. metropolitanus"
)

df <- df %>%
    mutate(
        genome = factor(
            genome,
            levels = genome_levels
        ),
        condition = factor(
            condition,
            levels = c(
                "control",
                "anhydrobiosis"
            )
        ),
        species = factor(
            species,
            levels = c(
                "experimentalis",
                "gadabouti"
            )
        )
    )

colors_condition <- c(
    control = "#1b9e77",
    anhydrobiosis = "#d95f02"
)

shapes_condition <- c(
    control = 16,
    anhydrobiosis = 17
)

condition_labels <- c(
    control = "Control",
    anhydrobiosis = "Anhydrobiosis"
)


# Means for each candidate genome


means <- df %>%
    group_by(species, genome) %>%
    summarise(
        mean_unique = mean(uniquely_mapped_pct, na.rm = TRUE),
        .groups = "drop"
    )



# Plot function


make_plot <- function(data, mean_data, species_name) {

    ggplot(
        data,
        aes(
            x = genome,
            y = uniquely_mapped_pct,
            color = condition,
            shape = condition
        )
    ) +

        # Individual RNA-seq libraries
        geom_point(
            position = position_jitter(
                width = 0.12,
                height = 0,
                seed = 42
            ),
            size = 3,
            alpha = 0.85
        ) +

        # Mean across libraries
        geom_point(
            data = mean_data,
            aes(
                x = genome,
                y = mean_unique
            ),
            inherit.aes = FALSE,
            shape = 95,
            size = 8,
            color = "black"
        ) +

        scale_color_manual(
            values = colors_condition,
            labels = c(
                control = "Control",
                anhydrobiosis = "Anhydrobiosis"
            )
        ) +

        scale_shape_manual(
            values = shapes_condition,
            labels = c(
                control = "Control",
                anhydrobiosis = "Anhydrobiosis"
            )
        ) +

        scale_x_discrete(
            labels = genome_labels
        ) +


        scale_y_continuous(
            limits = c(0, 23),
            breaks = seq(0, 25, 2.5),
            expand = expansion(mult = c(0, 0.03))
        )+

        labs(
            title = species_name,
            x = NULL,
            y = "Uniquely mapped reads (%)",
            color = "Condition",
            shape = "Condition"
        ) +

        theme_classic(base_size = 12) +

        theme(
            plot.title = element_text(
                face = "italic",
                hjust = 0.5,
                size = 12
            ),
            axis.text.x = element_text(
                angle = 45,
                hjust = 1
            ),
            legend.position = "bottom"
        )
}

# Separate species


exp_data <- df %>%
    filter(species == "experimentalis")

gad_data <- df %>%
    filter(species == "gadabouti")

exp_means <- means %>%
    filter(species == "experimentalis")

gad_means <- means %>%
    filter(species == "gadabouti")


# Generate panels


p_exp <- make_plot(
    exp_data,
    exp_means,
    "Pam. experimentalis"
)

p_gad <- make_plot(
    gad_data,
    gad_means,
    "Pam. gadabouti"
)



# Combine


combined_plot <- (
    p_exp + p_gad
) +
    plot_layout(
        guides = "collect"
    ) +
    plot_annotation(
        tag_levels = "A"
    ) &
    theme(
        legend.position = "bottom"
    )



# Save vector PDF


ggsave(
    filename = output_pdf,
    plot = combined_plot,
    width = 10,
    height = 5.5,
    units = "in",
    device = cairo_pdf
)

message("Saved figure to: ", output_pdf)