# workflow/rules/preprocessing/reference_assessment.smk

# Candidate reference genome assessment


if REFERENCE_ASSESSMENT_ENABLED:

    # 1. BUSCO completeness assessment


    if REFERENCE_ASSESSMENT.get("busco", {}).get("enabled", False):

        rule busco_candidate:
            input:
                fasta=lambda wc: CANDIDATE_CONFIG[wc.genome]["fasta"]
            output:
                result_dir=directory(f"{BUSCO_DIR}/busco_{{genome}}")
            params:
                lineage=REFERENCE_ASSESSMENT["busco"]["lineage"],
                mode=REFERENCE_ASSESSMENT["busco"].get("mode", "genome"),
                predictor=lambda wc: (
                    "--metaeuk"
                    if REFERENCE_ASSESSMENT["busco"].get("predictor") == "metaeuk"
                    else ""
                ),
                out=lambda wc: f"busco_{wc.genome}"
            threads: 8
            log:
                f"{LOGS}/reference_assessment/busco/{{genome}}.log"
            conda:
                "../../envs/preprocessing/busco.yaml"
            shell:
                """
                mkdir -p {BUSCO_DIR} $(dirname {log})

                busco \
                    -i {input.fasta:q} \
                    -o {params.out} \
                    -l {params.lineage:q} \
                    -m {params.mode:q} \
                    -c {threads} \
                    --out_path {BUSCO_DIR} \
                    {params.predictor} \
                    -f \
                    > {log} 2>&1
                """

        # 2. Consolidated BUSCO summary


        rule summarize_busco:
            input:
                expand(
                    f"{BUSCO_DIR}/busco_{{genome}}",
                    genome=CANDIDATES
                )
            output:
                summary=f"{REFERENCE_SUMMARY_DIR}/busco_summary.tsv"
            run:
                import glob
                import os

                os.makedirs(REFERENCE_SUMMARY_DIR, exist_ok=True)

                with open(output.summary, "w") as out:
                    out.write(
                        "genome\taccession\tcomplete_single\t"
                        "complete_duplicated\tfragmented\tmissing\t"
                        "n_buscos\tpredictor\n"
                    )

                    for genome in CANDIDATES:

                        accession = CANDIDATE_CONFIG[genome]["accession"]

                        result_dir = f"{BUSCO_DIR}/busco_{genome}"

                        summaries = glob.glob(
                            f"{result_dir}/**/short_summary*.txt",
                            recursive=True
                        )

                        if not summaries:
                            raise ValueError(
                                f"No BUSCO short summary found for {genome}"
                            )

                        summary_file = summaries[0]

                        values = {}

                        with open(summary_file) as handle:
                            for line in handle:
                                line = line.strip()

                                if "Complete and single-copy BUSCOs" in line:
                                    values["single"] = line.split()[0]

                                elif "Complete and duplicated BUSCOs" in line:
                                    values["duplicated"] = line.split()[0]

                                elif "Fragmented BUSCOs" in line:
                                    values["fragmented"] = line.split()[0]

                                elif "Missing BUSCOs" in line:
                                    values["missing"] = line.split()[0]

                                elif "Total BUSCO groups searched" in line:
                                    values["total"] = line.split()[0]

                        predictor = REFERENCE_ASSESSMENT["busco"].get(
                            "predictor",
                            "default"
                        )

                        out.write(
                            f"{genome}\t"
                            f"{accession}\t"
                            f"{values.get('single', 'NA')}\t"
                            f"{values.get('duplicated', 'NA')}\t"
                            f"{values.get('fragmented', 'NA')}\t"
                            f"{values.get('missing', 'NA')}\t"
                            f"{values.get('total', 'NA')}\t"
                            f"{predictor}\n"
                        )


    # 3. STAR indexes for candidate genomes

    if REFERENCE_ASSESSMENT.get("trial_mapping", {}).get("enabled", False):

        rule star_index_candidate:
            input:
                fasta=lambda wc: CANDIDATE_CONFIG[wc.genome]["fasta"]
            output:
                index=directory(f"{TRIAL_INDEX_DIR}/{{genome}}")
            threads: 8
            resources:
                mem_mb=160000,
                time_min=240
            log:
                f"{LOGS}/reference_assessment/star_index/{{genome}}.log"
            conda:
                "../../envs/alignment/star.yaml"
            shell:
                """
                mkdir -p {output.index} $(dirname {log})

                STAR \
                    --runMode genomeGenerate \
                    --genomeDir {output.index} \
                    --genomeFastaFiles {input.fasta:q} \
                    --runThreadN {threads} \
                    --limitGenomeGenerateRAM 160000000000 \
                    > {log} 2>&1
                """

        # 4. Trial mapping against each candidate genome
        rule star_trial_align:
            input:
                r1=f"{TRIMMED}/{{sample}}_R1_val_1.fq.gz",
                r2=f"{TRIMMED}/{{sample}}_R2_val_2.fq.gz",
                index=f"{TRIAL_INDEX_DIR}/{{genome}}"
            output:
                log_final=f"{TRIAL_ALIGN_DIR}/{{genome}}/{{sample}}_Log.final.out"
            params:
                opts=REFERENCE_ASSESSMENT["trial_mapping"].get(
                    "star_params",
                    "--outSAMtype None"
                ),
                prefix=lambda wc: (
                    f"{TRIAL_ALIGN_DIR}/{wc.genome}/tmp/{wc.sample}."
                )
            threads: 8
            resources:
                mem_mb=16000,
                time_min=240
            log:
                f"{LOGS}/reference_assessment/trial_mapping/{{genome}}/{{sample}}.log"
            conda:
                "../../envs/alignment/star.yaml"
            shell:
                """
                mkdir -p \
                    {TRIAL_ALIGN_DIR}/{wildcards.genome} \
                    {TRIAL_ALIGN_DIR}/{wildcards.genome}/tmp \
                    $(dirname {log})

                STAR \
                    --genomeDir {input.index:q} \
                    --readFilesIn {input.r1:q} {input.r2:q} \
                    --readFilesCommand zcat \
                    --runThreadN {threads} \
                    {params.opts} \
                    --outFileNamePrefix {params.prefix} \
                    > {log} 2>&1

                mv \
                    {params.prefix}Log.final.out \
                    {output.log_final}

                rm -f \
                    {params.prefix}Log.out \
                    {params.prefix}Log.progress.out \
                    {params.prefix}SJ.out.tab
                """

        # 5. Consolidated STAR trial mapping summary

        rule summarize_trial_mapping:
            input:
                expand(
                    f"{TRIAL_ALIGN_DIR}/{{genome}}/{{sample}}_Log.final.out",
                    genome=CANDIDATES,
                    sample=SAMPLES
                )
            output:
                summary=f"{REFERENCE_SUMMARY_DIR}/star_mapping_summary.tsv"
            params:
                trial_dir=TRIAL_ALIGN_DIR,
                genomes=CANDIDATES,
                samples=SAMPLES
            script:
                "../../scripts/preprocessing/summarize_trial_mapping.py"

        rule subsample_trial_reads:
            input:
                r1=f"{TRIMMED}/{{sample}}_R1_val_1.fq.gz",
                r2=f"{TRIMMED}/{{sample}}_R2_val_2.fq.gz"
            output:
                r1=f"{TRIAL_SUBSAMPLE_DIR}/{{sample}}_R1.subsampled.fq.gz",
                r2=f"{TRIAL_SUBSAMPLE_DIR}/{{sample}}_R2.subsampled.fq.gz"
            params:
                n=REFERENCE_ASSESSMENT["trial_mapping"].get(
                    "subsample_n",
                    2000000
                ),
                seed=REFERENCE_ASSESSMENT["trial_mapping"].get(
                    "subsample_seed",
                    42
                )
            threads: 2
            resources:
                mem_mb=4000,
                time_min=60
            log:
                f"{LOGS}/reference_assessment/subsampling/{{sample}}.log"
            conda:
                "../../envs/preprocessing/seqkit.yaml"
            shell:
                r"""
                mkdir -p {TRIAL_SUBSAMPLE_DIR} $(dirname {log})

                seqkit sample \
                    --number {params.n} \
                    --rand-seed {params.seed} \
                    {input.r1:q} \
                    -o {output.r1:q} \
                    2> {log}

                seqkit sample \
                    --number {params.n} \
                    --rand-seed {params.seed} \
                    {input.r2:q} \
                    -o {output.r2:q} \
                    2>> {log}

                echo "Subsampled {wildcards.sample}: {params.n} read pairs requested" >> {log}
                echo "Seed: {params.seed}" >> {log}
                echo -n "R1 reads: " >> {log}
                seqkit stats -T {output.r1:q} | tail -n 1 | cut -f4 >> {log}

                echo -n "R2 reads: " >> {log}
                seqkit stats -T {output.r2:q} | tail -n 1 | cut -f4 >> {log}
                """

        rule star_trial_align_subsampled:
            input:
                r1=f"{TRIAL_SUBSAMPLE_DIR}/{{sample}}_R1.subsampled.fq.gz",
                r2=f"{TRIAL_SUBSAMPLE_DIR}/{{sample}}_R2.subsampled.fq.gz",
                index=f"{TRIAL_INDEX_DIR}/{{genome}}"
            output:
                log_final=f"{TRIAL_ALIGN_SUBSAMPLED_DIR}/{{genome}}/{{sample}}_Log.final.out"
            params:
                opts=REFERENCE_ASSESSMENT["trial_mapping"].get(
                    "star_params",
                    "--outSAMtype None"
                ),
                prefix=lambda wc: (
                    f"{TRIAL_ALIGN_SUBSAMPLED_DIR}/"
                    f"{wc.genome}/tmp/{wc.sample}."
                )
            threads: 8
            resources:
                mem_mb=16000,
                time_min=120
            log:
                (
                    f"{LOGS}/reference_assessment/"
                    "trial_mapping_subsampled/{genome}/{sample}.log"
                )
            conda:
                "../../envs/alignment/star.yaml"
            shell:
                """
                mkdir -p \
                    {TRIAL_ALIGN_SUBSAMPLED_DIR}/{wildcards.genome} \
                    {TRIAL_ALIGN_SUBSAMPLED_DIR}/{wildcards.genome}/tmp \
                    $(dirname {log})

                STAR \
                    --genomeDir {input.index:q} \
                    --readFilesIn {input.r1:q} {input.r2:q} \
                    --readFilesCommand zcat \
                    --runThreadN {threads} \
                    {params.opts} \
                    --outFileNamePrefix {params.prefix} \
                    > {log} 2>&1

                mv \
                    {params.prefix}Log.final.out \
                    {output.log_final}

                rm -f \
                    {params.prefix}Log.out \
                    {params.prefix}Log.progress.out \
                    {params.prefix}SJ.out.tab
                """
                
        rule summarize_trial_mapping_subsampled:
            input:
                expand(
                    f"{TRIAL_ALIGN_SUBSAMPLED_DIR}/{{genome}}/{{sample}}_Log.final.out",
                    genome=CANDIDATES,
                    sample=SAMPLES
                )
            output:
                summary=(
                    f"{REFERENCE_SUMMARY_DIR}/"
                    "star_mapping_subsampled_summary.tsv"
                )
            params:
                trial_dir=TRIAL_ALIGN_SUBSAMPLED_DIR,
                genomes=CANDIDATES,
                samples=SAMPLES
            script:
                "../../scripts/preprocessing/summarize_trial_mapping.py"

        rule trial_mapping_subsampled_all:
            input:
                (
                    f"{REFERENCE_SUMMARY_DIR}/"
                    "star_mapping_subsampled_summary.tsv"
                )