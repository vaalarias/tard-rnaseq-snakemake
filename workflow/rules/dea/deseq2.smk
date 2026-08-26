rule save_raw_counts:
    input:
        samplesheet = config["samplesheet"],
        counts = expand(
            rules.featurecounts.output.counts,
            sample=SAMPLES
        )
    output:
        raw_counts = f"{DESEQ}/raw_counts_fc.rds"
    threads:
        1
    resources:
        mem_mb = 8000,
        time_min = 60
    conda:
        "../../envs/dea/r_dea.yaml"
    script:
        "../../scripts/dea/raw_counts_fc.R"

rule normalize_counts:
    input:
        raw_counts = f"{DESEQ}/raw_counts_fc.rds"
    output:
        vst  = f"{NORMALIZATION}/vst_matrix.rds",
        rlog = f"{NORMALIZATION}/rlog_matrix.rds",
        cpm  = f"{NORMALIZATION}/cpm_matrix.rds",
        tpm  = f"{NORMALIZATION}/tpm_matrix.rds"
    conda: 
        "../../envs/dea/r_dea.yaml"
    script:
        "../../scripts/normalization/normalization_counts.R"

rule run_deseq2:
    input:
        raw_counts = rules.save_raw_counts.output.raw_counts,
        config = "config/config.yaml"
    output:
        deseq_results = f"{DESEQ}/deseq_results_by_species.rds",
        exp_results = f"{DESEQ}/experimentalis_deseq2_results.csv",
        gad_results = f"{DESEQ}/gadabouti_deseq2_results.csv"
    params:
        outdir = DESEQ
    conda:
        "../../envs/dea/r_dea.yaml"
    script:
        "../../scripts/dea/deseq2.R"

rule summarize_degs:
    input:
        exp_results = rules.run_deseq2.output.exp_results,
        gad_results = rules.run_deseq2.output.gad_results,
        config = "config/config.yaml"
    output:
        species_summary = f"{DESEQ}/deg_summary_by_species.csv",
        shared_summary = f"{DESEQ}/shared_deg_summary.csv",
        gene_comparison = f"{DESEQ}/deg_gene_comparison.csv",
        gene_sets = f"{DESEQ}/deg_gene_sets.rds"
    resources:
        mem_mb = 4000,
        time_min = 30
    conda:
        "../../envs/dea/r_dea.yaml"
    script:
        "../../scripts/dea/summarize_degs.R"

rule annotate_genes:
    input:
        exp = rules.run_deseq2.output.exp_results,
        gad = rules.run_deseq2.output.gad_results,
        gaf = config["refs"]["gaf"]
    output:
        term2gene = f"{ANNOTATION}/TERM2GENE.rds",
        term2name = f"{ANNOTATION}/TERM2NAME.rds",
        annot_exp = f"{ANNOTATION}/experimentalis_annotated.rds",
        unan_exp  = f"{ANNOTATION}/experimentalis_unannotated.rds",
        annot_gad = f"{ANNOTATION}/gadabouti_annotated.rds",
        unan_gad  = f"{ANNOTATION}/gadabouti_unannotated.rds"
    conda: 
        "../../envs/dea/r_dea.yaml"
    script:
        "../../scripts/dea/annotate_genes.R"

rule go_enrichment:
    input:
        annot_exp = rules.annotate_genes.output.annot_exp,
        annot_gad = rules.annotate_genes.output.annot_gad,
        term2gene = rules.annotate_genes.output.term2gene,
        term2name = rules.annotate_genes.output.term2name,
        config = "config/config.yaml"
    output:
        exp_up_pdf = f"{PLOTS}/enrichment/experimentalis_up_go.pdf",
        exp_down_pdf = f"{PLOTS}/enrichment/experimentalis_down_go.pdf",
        gad_up_pdf = f"{PLOTS}/enrichment/gadabouti_up_go.pdf",
        gad_down_pdf = f"{PLOTS}/enrichment/gadabouti_down_go.pdf"
    resources:
        mem_mb = 8000,
        time_min = 60
    conda:
        "../../envs/dea/r_dea.yaml"
    script:
        "../../scripts/dea/go_enrichment.R"

rule summarize_annotations:
    input:
        annot_exp = rules.annotate_genes.output.annot_exp,
        unan_exp = rules.annotate_genes.output.unan_exp,
        annot_gad = rules.annotate_genes.output.annot_gad,
        unan_gad = rules.annotate_genes.output.unan_gad,
        config = "config/config.yaml"
    output:
        annotation_summary = f"{DESEQ}/deg_annotation_summary.csv",
        shared_summary = f"{DESEQ}/shared_deg_annotation_summary.csv",
        shared_details = f"{DESEQ}/shared_deg_annotation_details.csv"
    conda:
        "../../envs/dea/r_dea.yaml"
    script:
        "../../scripts/dea/summarize_annotations.R"

rule qc_summary:
    input:
        samplesheet = config["samplesheet"],
        raw_counts = f"{DESEQ}/raw_counts_fc.rds",
        config = "config/config.yaml",
        align_dir = f"{ALIGN}/logs"
    output:
        "results/plots/qc_summary.pdf"
    conda:
        "../../envs/dea/r_dea.yaml"
    script:
        "../../scripts/plots/qc_summary.R"
