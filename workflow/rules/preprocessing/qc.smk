rule fastqc_raw:
    input:
        fq=f"{RAW_FASTQ}/{{sample}}_R{{read}}.fastq.gz"
    output:
        html=f"{FASTQC_RAW}/{{sample}}_R{{read}}_fastqc.html",
        zip=f"{FASTQC_RAW}/{{sample}}_R{{read}}_fastqc.zip"
    threads: 2
    log:
        str(LOGS / "fastqc_raw" / "{sample}_R{read}.log")
    conda:
        "../../envs/preprocessing/fastqc.yaml"
    shell:
        """
        mkdir -p {FASTQC_RAW} $(dirname {log})
        fastqc -t {threads} -o {FASTQC_RAW} {input.fq:q} > {log} 2>&1
        """


rule trim_galore:
    input:
        r1=f"{RAW_FASTQ}/{{sample}}_R1.fastq.gz",
        r2=f"{RAW_FASTQ}/{{sample}}_R2.fastq.gz"
    output:
        r1=f"{TRIMMED}/{{sample}}_R1_val_1.fq.gz",
        r2=f"{TRIMMED}/{{sample}}_R2_val_2.fq.gz",
        report_r1=f"{TRIMMED}/{{sample}}_R1.fastq.gz_trimming_report.txt",
        report_r2=f"{TRIMMED}/{{sample}}_R2.fastq.gz_trimming_report.txt"
    threads: 2
    log:
        str(LOGS / "trimgalore" / "{sample}.log")
    conda:
        "../../envs/preprocessing/trimgalore.yaml"
    shell:
        """
        mkdir -p {TRIMMED} $(dirname {log})

        trim_galore --paired --cores {threads} \
            --output_dir {TRIMMED} \
            {input.r1:q} {input.r2:q} > {log} 2>&1
        """


rule fastqc_trimmed:
    input:
        fq=f"{TRIMMED}/{{sample}}_R{{read}}_val_{{read}}.fq.gz"
    output:
        html=f"{FASTQC_TRIMMED}/{{sample}}_R{{read}}_val_{{read}}_fastqc.html",
        zip=f"{FASTQC_TRIMMED}/{{sample}}_R{{read}}_val_{{read}}_fastqc.zip"
    threads: 2
    log:
        str(LOGS / "fastqc_trimmed" / "{sample}_R{read}.log")
    conda:
        "../../envs/preprocessing/fastqc.yaml"
    shell:
        """
        mkdir -p {FASTQC_TRIMMED} $(dirname {log})
        fastqc -t {threads} -o {FASTQC_TRIMMED} {input.fq:q} > {log} 2>&1
        """