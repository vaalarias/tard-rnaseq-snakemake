rule summarize_star_unmapped_reads:
    input:
        r1 = (
            f"{ALIGN}/star/unmapped/"
            "{sample}_R1_unmapped.fq"
        ),
        r2 = (
            f"{ALIGN}/star/unmapped/"
            "{sample}_R2_unmapped.fq"
        )
    output:
        summary = (
            f"{BLAST}/unmapped_read_sequences/"
            f"{{sample}}.unmapped_sequence_summary.tsv"
        ),
        fasta = (
            f"{BLAST}/unmapped_read_sequences/"
            f"{{sample}}.unmapped_sequences.fa"
        )
    params:
        head_n = config["unmapped_reads"]["head_n"],
        fraction = config["unmapped_reads"]["subsample_fraction"],
        min_length = config["unmapped_reads"]["min_length"],
        max_n_fraction = config["unmapped_reads"]["max_n_fraction"],
        seed = config["unmapped_reads"]["random_seed"],
        max_base_fraction = (config["unmapped_reads"]["max_base_fraction"]),
    log:
        (
            f"{LOGS}/blast/unmapped_reads/"
            f"{{sample}}.log"
        )
    conda:
        "../../envs/contamination_screening/py_script.yml"
    script:
        "../../scripts/contamination_screening/summarize_unmapped_reads.py"

BLAST_SAMPLES = config["blast"]["samples"]


rule combine_unmapped_blast_queries:
    input:
        [
            (
                f"{BLAST}/unmapped_read_sequences/"
                f"{sample}.unmapped_sequences.fa"
            )
            for sample in BLAST_SAMPLES
        ]
    output:
        fasta = (
            f"{BLAST}/unmapped_read_sequences/"
            "all_samples.unmapped_sequences.fa"
        )
    shell:
        r"""
        mkdir -p $(dirname {output.fasta:q})
        cat {input:q} > {output.fasta:q}
        test -s {output.fasta:q}
        """


rule blast_unmapped_sequences:
    input:
        fasta = (
            rules.combine_unmapped_blast_queries
            .output.fasta
        )
    output:
        blast = (
            f"{BLAST}/unmapped_read_sequences/"
            "all_samples.blast_nt.tsv"
        )
    params:
        database = config["blast"]["database"],
        task = config["blast"]["task"],
        evalue = config["blast"]["evalue"],
        max_target_seqs = (
            config["blast"]["max_target_seqs"]
        ),
        chunk_size = (
            config["blast"].get(
                "chunk_size",
                10
            )
        ),
        max_retries = (
            config["blast"].get(
                "max_retries",
                3
            )
        ),
        retry_wait_seconds = (
            config["blast"].get(
                "retry_wait_seconds",
                60
            )
        ),
        attempt_timeout_seconds = (
            config["blast"].get(
                "attempt_timeout_seconds",
                1800
            )
        )
    resources:
        mem_mb = 4000,
        time_min = 720,
        http = 1
    conda:
        "../../envs/contamination_screening/blast.yaml"
    log:
        f"{LOGS}/blast/unmapped_reads/blast_nt.log"
    script:
        "../../scripts/contamination_screening/blast_unmapped_sequences.py"

rule summarize_unmapped_blast:
    input:
        blast = rules.blast_unmapped_sequences.output.blast,
        fasta = rules.combine_unmapped_blast_queries.output.fasta
    output:
        best_hits = (
            f"{BLAST}/unmapped_read_sequences/"
            "unmapped_blast_best_hits.tsv"
        ),
        taxon_summary = (
            f"{BLAST}/unmapped_read_sequences/"
            "unmapped_blast_taxon_summary.tsv"
        )
    conda:
        "../../envs/contamination_screening/py_script.yml"
    script:
        "../../scripts/contamination_screening/summarize_unmapped_blast.py"
