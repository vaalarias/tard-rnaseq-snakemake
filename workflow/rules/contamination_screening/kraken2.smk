rule kraken2_unmapped:
    input:
        r1=f"{ALIGN}/star/unmapped/{{sample}}_R1_unmapped.fq",
        r2=f"{ALIGN}/star/unmapped/{{sample}}_R2_unmapped.fq"
    output:
        report=f"{KRAKEN}/{{sample}}_unmapped.kraken2.report",
        classified_r1=f"{KRAKEN}/{{sample}}_unmapped.classified_1.fastq",
        classified_r2=f"{KRAKEN}/{{sample}}_unmapped.classified_2.fastq",
        unclassified_r1=f"{KRAKEN}/{{sample}}_unmapped.unclassified_1.fastq",
        unclassified_r2=f"{KRAKEN}/{{sample}}_unmapped.unclassified_2.fastq"
    params:
        db=config["kraken2"]["db"],
        confidence=config["kraken2"]["confidence"],
        minimum_hit_group=config["kraken2"].get("minimum_hit_group", 3)
    log:
        f"{LOGS}/kraken2/{{sample}}_unmapped.kraken2.log"
    conda:
        "../../envs/contamination_screening/kraken.yaml"
    threads: 12
    shell:
        """
        mkdir -p {KRAKEN} $(dirname {log})

        kraken2 \
            --db {params.db} \
            --paired \
            --report {output.report} \
            --confidence {params.confidence} \
            --classified-out {KRAKEN}/{wildcards.sample}_unmapped.classified#.fastq \
            --unclassified-out {KRAKEN}/{wildcards.sample}_unmapped.unclassified#.fastq \
            --minimum-hit-group {params.minimum_hit_group} \
            --threads {threads} \
            {input.r1:q} {input.r2:q} > {log} 2>&1
        """

rule kraken2_trimmed:
    input:
        r1=f"{TRIMMED}/{{sample}}_R1_val_1.fq.gz",
        r2=f"{TRIMMED}/{{sample}}_R2_val_2.fq.gz"
    output:
        report=f"{KRAKEN}/{{sample}}_trimmed.kraken2.report",
        classified_r1=f"{KRAKEN}/{{sample}}_trimmed.classified_1.fastq",
        classified_r2=f"{KRAKEN}/{{sample}}_trimmed.classified_2.fastq",
        unclassified_r1=f"{KRAKEN}/{{sample}}_trimmed.unclassified_1.fastq",
        unclassified_r2=f"{KRAKEN}/{{sample}}_trimmed.unclassified_2.fastq"
    params:
        db=config["kraken2"]["db"],
        confidence=config["kraken2"]["confidence"],
        minimum_hit_group=config["kraken2"].get("minimum_hit_group", 3)
    log:
        f"{LOGS}/kraken2/{{sample}}_trimmed.kraken2.log"
    threads: 12
    conda:
        "../../envs/contamination_screening/kraken.yaml"
    shell:
        """
        mkdir -p {KRAKEN} $(dirname {log})

        kraken2 \
            --gzip-compressed \
            --db {params.db} \
            --paired \
            --report {output.report} \
            --confidence {params.confidence} \
            --classified-out {KRAKEN}/{wildcards.sample}_trimmed.classified#.fastq \
            --unclassified-out {KRAKEN}/{wildcards.sample}_trimmed.unclassified#.fastq \
            --minimum-hit-group {params.minimum_hit_group} \
            --threads {threads} \
            {input.r1:q} {input.r2:q} > {log} 2>&1
        """

rule bracken_unmapped:
    input:
        report=f"{KRAKEN}/{{sample}}_unmapped.kraken2.report"
    output:
        bracken=f"{KRAKEN}/{{sample}}_unmapped.bracken"
    params:
        level=config["kraken2"].get("bracken_level", "S"),
        readlen=config["kraken2"].get("bracken_readlen", 100),
        db=config["kraken2"]["db"]
    log:
        f"{LOGS}/kraken2/{{sample}}_unmapped.bracken.log"
    conda:
        "../../envs/contamination_screening/kraken.yaml"
    shell:
        """
        mkdir -p {KRAKEN} $(dirname {log})

        bracken \
            -d {params.db} \
            -i {input.report} \
            -o {output.bracken} \
            -r {params.readlen} \
            -l {params.level} > {log} 2>&1
        """

rule bracken_trimmed:
    input:
        report=f"{KRAKEN}/{{sample}}_trimmed.kraken2.report"
    output:
        bracken=f"{KRAKEN}/{{sample}}_trimmed.bracken"
    params:
        db=config["kraken2"]["db"],
        level=config["kraken2"].get("bracken_level", "S"),
        readlen=config["kraken2"].get("bracken_readlen", 100)
    log:
        f"{LOGS}/kraken2/{{sample}}_trimmed.bracken.log"
    conda:
        "../../envs/contamination_screening/kraken.yaml"
    threads: 4
    shell:
        """
        mkdir -p {KRAKEN} $(dirname {log})

        bracken \
            -d {params.db} \
            -i {input.report} \
            -o {output.bracken} \
            -r {params.readlen} \
            -l {params.level} > {log} 2>&1
        """