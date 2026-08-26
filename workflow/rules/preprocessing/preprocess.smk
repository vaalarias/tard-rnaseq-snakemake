rule merge_replicates:
    input:
        unpack(input_merge)
    output:
        r1=f"{RAW_FASTQ}/{{sample}}_R1.fastq.gz",
        r2=f"{RAW_FASTQ}/{{sample}}_R2.fastq.gz"
    log:
        str(LOGS / "merge_replicates" / "{sample}.log")
    shell:
        """
        mkdir -p {RAW_FASTQ} $(dirname {log})

        zcat {input.r1:q} | gzip -c > {output.r1}
        zcat {input.r2:q} | gzip -c > {output.r2}
        """