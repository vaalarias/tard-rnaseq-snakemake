# workflow/rules/targets.smk
# workflow/rules/targets.smk


def targets():

    ans = []

    # ========================================================
    # Main analysis outputs
    # ========================================================

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

    # ========================================================
    # BLAST analysis of unmapped reads
    # ========================================================

    ans.extend([
        f"{BLAST}/unmapped_read_sequences/unmapped_blast_best_hits.tsv",
        f"{BLAST}/unmapped_read_sequences/unmapped_blast_taxon_summary.tsv",
    ])

    # ========================================================
    # Candidate reference genome assessment
    # ========================================================

    if REFERENCE_ASSESSMENT_ENABLED:

        # BUSCO
        if REFERENCE_ASSESSMENT.get("busco", {}).get("enabled", False):

            ans.extend(
                expand(
                    f"{BUSCO_DIR}/busco_{{genome}}",
                    genome=CANDIDATES
                )
            )

            ans.append(
                f"{REFERENCE_SUMMARY_DIR}/busco_summary.tsv"
            )

        # STAR trial mapping - subsampled reads
        if REFERENCE_ASSESSMENT.get(
            "trial_mapping", {}
        ).get("enabled", False):

            ans.extend(
                expand(
                    f"{TRIAL_ALIGN_SUBSAMPLED_DIR}/{{genome}}/{{sample}}_Log.final.out",
                    genome=CANDIDATES,
                    sample=SAMPLES
                )
            )

            ans.append(
                f"{REFERENCE_SUMMARY_DIR}/star_mapping_subsampled_summary.tsv"
            )



    return ans