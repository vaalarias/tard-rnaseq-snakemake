PLOT_SPECIES = config["plots"]["species"]
EPI_FAMILIES = config["plots"]["epigenetic_families"]

rule figures:
    input:
        plot_targets()

rule build_species_de_objects:
    input:
        raw_counts = f"{DESEQ}/raw_counts_fc.rds"
    output:
        dds_list = f"{DESEQ}/dds_list.rds",
        vst_list = f"{DESEQ}/vst_list.rds"
    resources:
        mem_mb = 16000,
        time_min = 120
    conda:
        "../../envs/plots/plots.yml"
    script:
        "../../scripts/plots/build_species_objects.R"

rule extract_histone_core:
    input:
        annotation = GTF
    output:
        histones = f"{DESEQ}/histone_core.tsv"
    shell:
        r"""
        mkdir -p $(dirname {output.histones:q})

        awk -F '\t' '
        BEGIN {{
            OFS = "\t"
        }}

        $0 !~ /^#/ &&
        tolower($9) ~ /product ".*histone (h1|h2|h3|h4)/ {{
            gene = ""
            product = ""
            protein = ""

            if (match($9, /gene_id "([^"]+)"/, a))
                gene = a[1]

            if (match($9, /product "([^"]+)"/, b))
                product = b[1]

            if (match($9, /protein_id "([^"]+)"/, c))
                protein = c[1]

            if (gene != "" && product != "" && !seen[gene]++)
                print gene, product, protein
        }}
        ' {input.annotation:q} > {output.histones:q}

        test -s {output.histones:q}
        """

rule extract_epigenetic_elements:
    input:
        annotation = GTF
    output:
        regulators = (
            f"{OUTDIR}/annotations/"
            "epigenetic_regulators.tsv"
        ),
        epitranscriptomic = (
            f"{OUTDIR}/annotations/"
            "epitranscriptomic.txt"
        )
    shell:
        r"""
        mkdir -p $(dirname {output.regulators:q})

        # Chromatin and histone regulators
        awk -F '\t' '
        BEGIN {{
            OFS = "\t"
        }}

        $0 !~ /^#/ &&
        tolower($9) ~ /(polycomb|prc1|prc2|ezh1|ezh2|suz12|eed|ring1|bmi1|histone.*demethylase|jumonji|kdm[0-9]|histone.*methyltransferase|set domain|kmt[0-9]|smyd|histone.*deacetylase|hdac[0-9]|sirtuin|histone.*acetyltransferase|gnat|myst)/ {{
            gene = ""
            product = ""
            protein = ""

            if (match($9, /gene_id "([^"]+)"/, a))
                gene = a[1]

            if (match($9, /product "([^"]+)"/, b))
                product = b[1]

            if (match($9, /protein_id "([^"]+)"/, c))
                protein = c[1]

            if (gene != "" && product != "" && !seen[gene]++)
                print gene, product, protein
        }}
        ' {input.annotation:q} > {output.regulators:q}

        # RNA-modification and epitranscriptomic elements
        awk -F '\t' '
        BEGIN {{
            OFS = "\t"
        }}

        $0 !~ /^#/ &&
        tolower($9) ~ /(rna.*methyltransferase|trna.*methyltransferase|rrna.*methyltransferase|pseudouridine synthase|mettl3|mettl14|wtap|virma|kiaa1429|rbm15|zc3h13|fto|alkbh5|ythdf|ythdc|igf2bp|nsun|trdmt1|dnmt2|adar|apobec)/ {{
            gene = ""
            product = ""
            protein = ""

            if (match($9, /gene_id "([^"]+)"/, a))
                gene = a[1]

            if (match($9, /product "([^"]+)"/, b))
                product = b[1]

            if (match($9, /protein_id "([^"]+)"/, c))
                protein = c[1]

            if (gene != "" && product != "" && !seen[gene]++)
                print gene, product, protein
        }}
        ' {input.annotation:q} > {output.epitranscriptomic:q}

        test -s {output.regulators:q}
        test -s {output.epitranscriptomic:q}
        """

rule build_qc_summary:
    input:
        raw_counts = f"{DESEQ}/raw_counts_fc.rds",
        plot_metadata = "config/plot_metadata.yaml"
    output:
        qc_summary = f"{OUTDIR}/qc_summary.csv"
    conda:
        "../../envs/plots/plots.yml"
    script:
        "../../scripts/plots/build_qc_summary.R"


rule plot_qc_combined_boxplots:
    input:
        qc_summary = rules.build_qc_summary.output.qc_summary
    output:
        plot = f"{PLOTS}/qc_combined_boxplots.pdf"
    conda:
        "../../envs/plots/plots.yml"
    script:
        "../../scripts/plots/plot_qc_summary.R"


rule plot_pca:
    input:
        vst = f"{NORMALIZATION}/vst_matrix.rds",
        vst_list = rules.build_species_de_objects.output.vst_list,
        raw_counts = f"{DESEQ}/raw_counts_fc.rds",
        plot_metadata = "config/plot_metadata.yaml"
    output:
        total = f"{PLOTS}/pca_vst.pdf",
        species = expand(f"{PLOTS}/pca_{{species}}.pdf", species=PLOT_SPECIES)
    conda:
        "../../envs/plots/plots.yml"
    script:
        "../../scripts/plots/plot_pca.R"


rule plot_sample_correlation:
    input:
        vst = f"{NORMALIZATION}/vst_matrix.rds",
        raw_counts = f"{DESEQ}/raw_counts_fc.rds"
    output:
        plot = f"{PLOTS}/heatmap_sample_correlation.pdf"
    conda:
        "../../envs/plots/plots.yml"
    script:
        "../../scripts/plots/plot_sample_correlation.R"


rule plot_deg_heatmaps:
    input:
        deseq_results = f"{DESEQ}/deseq_results_by_species.rds",
        vst_list = rules.build_species_de_objects.output.vst_list,
        raw_counts = f"{DESEQ}/raw_counts_fc.rds"
    output:
        expand(f"{PLOTS}/heatmap_DEGs_{{species}}.pdf", species=PLOT_SPECIES)
    resources:
        mem_mb = 8000,
        time_min = 60
    conda:
        "../../envs/plots/plots.yml"
    script:
        "../../scripts/plots/plot_deg_heatmaps.R"


rule plot_venn_degs:
    input:
        deseq_results = f"{DESEQ}/deseq_results_by_species.rds"
    output:
        plot = f"{PLOTS}/venn_DEGs.pdf"
    conda:
        "../../envs/plots/plots.yml"
    script:
        "../../scripts/plots/plot_venn_degs.R"


rule plot_volcanoes:
    input:
        deseq_results = f"{DESEQ}/deseq_results_by_species.rds"
    output:
        expand(f"{PLOTS}/volcano_{{species}}.pdf", species=["Experimentalis", "Gadabouti"])
    conda:
        "../../envs/plots/plots.yml"
    script:
        "../../scripts/plots/plot_volcanoes.R"


rule plot_histone_figures:
    input:
        tpm = f"{NORMALIZATION}/tpm_matrix.rds",
        raw_counts = f"{DESEQ}/raw_counts_fc.rds",
        histones = f"{DESEQ}/histone_core.tsv"
    output:
        z_uncollapsed = expand(
            f"{PLOTS}/heatmap_histones_Zscores_{{species}}no_collapse.pdf",
            species=PLOT_SPECIES
        ),
        tpm_species = expand(
            f"{PLOTS}/heatmap_histones_TPM_{{species}}.pdf",
            species=PLOT_SPECIES
        ),
        logtpm_species = expand(
            f"{PLOTS}/heatmap_histones_logTPM_{{species}}.pdf",
            species=PLOT_SPECIES
        ),
        z_merged = expand(
            f"{PLOTS}/heatmap_histones_Zscores_{{species}}merged.pdf",
            species=PLOT_SPECIES
        ),
        polya_no = expand(
            f"{PLOTS}/heatmap_histones_no_Zscores_{{species}}.pdf",
            species=PLOT_SPECIES
        ),
        polya_yes = expand(
            f"{PLOTS}/heatmap_histones_yes_Zscores_{{species}}.pdf",
            species=PLOT_SPECIES
        ),
        tpm_all = f"{PLOTS}/heatmap_histones_TPM.pdf",
        z_all = f"{PLOTS}/heatmap_histones_Zscores.pdf"
    resources:
        mem_mb = 8000,
        time_min = 120
    conda:
        "../../envs/plots/plots.yml"
    script:
        "../../scripts/plots/plot_histones.R"

rule run_dea_method_comparison:
    input:
        raw_counts = f"{DESEQ}/raw_counts_fc.rds",
        config = "config/config.yaml"
    output:
        deseq2 = f"{DESEQ}/methods/deseq2_results.rds",
        edger = f"{DESEQ}/methods/edger_results.rds",
        limma = f"{DESEQ}/methods/limma_voom_results.rds",
        diagnostics = f"{DESEQ}/methods/method_diagnostics.rds",
        summary = f"{DESEQ}/methods/method_overlap_summary.csv"
    resources:
        mem_mb = 16000,
        time_min = 180
    conda:
        "../../envs/plots/plots.yml"
    script:
        "../../scripts/plots/run_dea_methods.R"


rule plot_dea_method_comparison:
    input:
        deseq2 = rules.run_dea_method_comparison.output.deseq2,
        edger = rules.run_dea_method_comparison.output.edger,
        limma = rules.run_dea_method_comparison.output.limma,
        diagnostics = rules.run_dea_method_comparison.output.diagnostics,
        config = "config/config.yaml"
    output:
        venn = expand(f"{PLOTS}/dea_methods/venn_methods_{{species}}.pdf", species=PLOT_SPECIES),
        density = expand(f"{PLOTS}/dea_methods/log2fc_density_{{species}}.pdf", species=PLOT_SPECIES),
        bcv = expand(f"{PLOTS}/dea_methods/edgeR_BCV_{{species}}.pdf", species=PLOT_SPECIES),
        voom = expand(f"{PLOTS}/dea_methods/limma_voom_{{species}}.pdf", species=PLOT_SPECIES),
        ma = expand(
            f"{PLOTS}/dea_methods/MA_{{method}}_{{species}}.pdf",
            method=["DESeq2", "edgeR", "limma_voom"],
            species=PLOT_SPECIES
        )
    resources:
        mem_mb = 8000,
        time_min = 120
    conda:
        "../../envs/plots/plots.yml"
    script:
        "../../scripts/plots/plot_dea_methods.R"


rule plot_cross_species_scatter:
    input:
        vst_list = rules.build_species_de_objects.output.vst_list,
        tpm = f"{NORMALIZATION}/tpm_matrix.rds",
        raw_counts = f"{DESEQ}/raw_counts_fc.rds"
    output:
        control_vst = f"{PLOTS}/scatter_control_vst.pdf",
        control_tpm = f"{PLOTS}/scatter_control_tpm.pdf",
        anh_vst = f"{PLOTS}/scatter_anhydrobiosis_vst.pdf",
        anh_tpm = f"{PLOTS}/scatter_anhydrobiosis_tpm.pdf"
    conda:
        "../../envs/plots/plots.yml"
    script:
        "../../scripts/plots/plot_cross_species_scatter.R"


rule plot_normalization_qc:
    input:
        raw_counts = f"{DESEQ}/raw_counts_fc.rds",
        tpm = f"{NORMALIZATION}/tpm_matrix.rds",
        vst = f"{NORMALIZATION}/vst_matrix.rds"
    output:
        tpm_hist = f"{PLOTS}/normalization/tpm_distribution_histogram.pdf",
        tpm_box = f"{PLOTS}/normalization/tpm_distribution_by_sample.pdf",
        vst_hist = f"{PLOTS}/normalization/vst_distribution_histogram.pdf",
        vst_box = f"{PLOTS}/normalization/vst_distribution_by_sample.pdf"
    conda:
        "../../envs/plots/plots.yml"
    script:
        "../../scripts/plots/plot_normalization_qc.R"


rule plot_extended_feature_heatmaps:
    input:
        raw_counts = f"{DESEQ}/raw_counts_fc.rds",
        tpm = f"{NORMALIZATION}/tpm_matrix.rds",
        vst = f"{NORMALIZATION}/vst_matrix.rds",
        deseq_results = f"{DESEQ}/deseq_results_by_species.rds",
        histones = f"{DESEQ}/histone_core.tsv",
        epigenetic_regulators = (
            f"{OUTDIR}/annotations/epigenetic_regulators.tsv"
        ),
        epitranscriptomic = (
            f"{OUTDIR}/annotations/epitranscriptomic.txt"
        )
    output:
        histone_dir = directory(f"{PLOTS}/heatmaps/histones"),
        epigenetic_dir = directory(f"{PLOTS}/heatmaps/epigenetic_elements")
    resources:
        mem_mb = 16000,
        time_min = 240
    conda:
        "../../envs/plots/plots.yml"
    script:
        "../../scripts/plots/plot_extended_feature_heatmaps.R"


rule plot_go_interest:
    input:
        raw_counts = f"{DESEQ}/raw_counts_fc.rds",
        tpm = f"{NORMALIZATION}/tpm_matrix.rds",
        deseq_results = f"{DESEQ}/deseq_results_by_species.rds",
        annot_exp = f"{ANNOTATION}/experimentalis_annotated.rds",
        annot_gad = f"{ANNOTATION}/gadabouti_annotated.rds",
        term2gene = f"{ANNOTATION}/TERM2GENE.rds",
        term2name = f"{ANNOTATION}/TERM2NAME.rds",
        gaf = config["refs"]["gaf"]
    output:
        go_dir = directory(f"{PLOTS}/heatmaps/go_interest")
    resources:
        mem_mb = 16000,
        time_min = 240
    conda:
        "../../envs/plots/plots.yml"
    script:
        "../../scripts/plots/plot_go_interest.R"

rule plot_quantification_method_comparison:
    input:
        samplesheet = "data/samplesheet.csv",
        gtf = (
            "config/ref/"
            "GCF_019649055.1_Prichtersi_v1.0_genomic.gtf"
        ),
        salmon = [
            f"{ALIGN}/salmon/{sample}/quant.sf"
            for sample in SAMPLES
        ],
        featurecounts = [
            f"{COUNTS}/{sample}_counts.txt"
            for sample in SAMPLES
        ]
    output:
        method_plot = (
            f"{PLOTS}/method_comparison/"
            "figure_03b_method_comparison.pdf"
        ),
        overlap_plot = (
            f"{PLOTS}/method_comparison/"
            "figure_03a_feature_overlap.pdf"
        ),
        combined_plot = (
            f"{PLOTS}/method_comparison/"
            "figure_03_quantification_comparison.pdf"
        ),
        sample_table = (
            f"{PLOTS}/method_comparison/"
            "detected_features_by_sample.tsv"
        ),
        overlap_table = (
            f"{PLOTS}/method_comparison/"
            "feature_overlap_by_group.tsv"
        )
    params:
        samples = SAMPLES
    threads:
        1
    resources:
        mem_mb = 8000,
        time_min = 240
    conda:
        "../../envs/dea/r_dea.yaml"
    log:
        f"{LOGS}/plots/quantification_method_comparison.log"
    script:
        "../../scripts/plots/plot_quantification_method_comparison.R"