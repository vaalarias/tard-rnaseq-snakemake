def targets():

    ans = []

    ans.extend([
        "results/multiqc/multiqc_report.html",

        "results/normalization/vst_matrix.rds",
        "results/normalization/rlog_matrix.rds",
        "results/normalization/cpm_matrix.rds",
        "results/normalization/tpm_matrix.rds",

        "results/deseq/deseq_results_by_species.rds",

        "results/plots/enrichment/experimentalis_up_go.pdf",
        "results/plots/enrichment/experimentalis_down_go.pdf",
        "results/plots/enrichment/gadabouti_up_go.pdf",
        "results/plots/enrichment/gadabouti_down_go.pdf",
    ])

    ans.extend([
    f"{BLAST}/unmapped_read_sequences/"
    "unmapped_blast_best_hits.tsv",

    f"{BLAST}/unmapped_read_sequences/"
    "unmapped_blast_taxon_summary.tsv",
    ])

    ans.extend(plot_targets())

    return ans

