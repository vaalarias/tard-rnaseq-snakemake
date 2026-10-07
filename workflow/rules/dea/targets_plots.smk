def plot_targets():

    return [
        "results/qc_summary.csv",
        "results/plots/qc_combined_boxplots.pdf",
        "results/plots/pca_panel.pdf",
        "results/plots/heatmap_sample_correlation.pdf",
        # Reference-genome assessment
        "results/plots/reference_assessment/"
        "reference_trial_mapping.pdf",


        # Individual DEG plots
        "results/plots/heatmap_DEGs_experimentalis.pdf",
        "results/plots/heatmap_DEGs_gadabouti.pdf",
        "results/plots/volcano_Experimentalis.pdf",
        "results/plots/volcano_Gadabouti.pdf",

        # Combined volcano and heatmap panels
        "results/plots/DEG_panel_experimentalis.pdf",
        "results/plots/DEG_panel_gadabouti.pdf",

        # Shared and species-specific DEG analysis
        "results/deseq/shared_response/de_gene_classification.csv",
        "results/deseq/shared_response/de_response_summary.csv",
        "results/deseq/shared_response/de_overlap_summary.csv",
        "results/plots/shared_and_species_specific_DEGs.pdf",

        # Previous DEG overlap plot
        "results/plots/venn_DEGs.pdf",

        *expand(
            "results/plots/dea_methods/venn_methods_{species}.pdf",
            species=[
                "experimentalis",
                "gadabouti"
            ]
        ),

        *expand(
            "results/plots/dea_methods/log2fc_density_{species}.pdf",
            species=[
                "experimentalis",
                "gadabouti"
            ]
        ),

        *expand(
            "results/plots/dea_methods/edgeR_BCV_{species}.pdf",
            species=[
                "experimentalis",
                "gadabouti"
            ]
        ),

        *expand(
            "results/plots/dea_methods/limma_voom_{species}.pdf",
            species=[
                "experimentalis",
                "gadabouti"
            ]
        ),

        *expand(
            "results/plots/dea_methods/MA_{method}_{species}.pdf",
            method=[
                "DESeq2",
                "edgeR",
                "limma_voom"
            ],
            species=[
                "experimentalis",
                "gadabouti"
            ]
        ),

        "results/plots/scatter_control_vst.pdf",
        "results/plots/scatter_control_tpm.pdf",
        "results/plots/scatter_anhydrobiosis_vst.pdf",
        "results/plots/scatter_anhydrobiosis_tpm.pdf",

        "results/plots/normalization/"
        "tpm_distribution_histogram.pdf",

        "results/plots/normalization/"
        "tpm_distribution_by_sample.pdf",

        "results/plots/normalization/"
        "vst_distribution_histogram.pdf",

        "results/plots/normalization/"
        "vst_distribution_by_sample.pdf",

        "results/plots/heatmaps/histones",
        "results/plots/heatmaps/epigenetic_elements",
        "results/plots/heatmaps/go_interest",

        "results/plots/method_comparison/"
        "figure_03a_feature_overlap.pdf",

        "results/plots/method_comparison/"
        "figure_03b_method_comparison.pdf",

        "results/plots/method_comparison/"
        "figure_03_quantification_comparison.pdf",

        "results/plots/method_comparison/"
        "detected_features_by_sample.tsv",

        "results/plots/method_comparison/"
        "feature_overlap_by_group.tsv",

        
    ]