rule multiqc:
    input:
        expand(f"{FASTQC_RAW}/{{sample}}_R{{read}}_fastqc.html", sample=SAMPLES, read=[1,2]),
        expand(f"{FASTQC_TRIMMED}/{{sample}}_R{{read}}_val_{{read}}_fastqc.html", read=[1,2], sample=SAMPLES),
        expand(f"{ALIGN}/star/{{sample}}.bam", sample=SAMPLES),
        expand(f"{COUNTS}/{{sample}}_counts.txt", sample=SAMPLES),
        expand(f"{ALIGN}/salmon/{{sample}}", sample=SAMPLES),
        expand(f"{KRAKEN}/{{sample}}_trimmed.bracken", sample=SAMPLES),
        expand(f"{KRAKEN}/{{sample}}_unmapped.bracken", sample=SAMPLES)

    output:
        directory(f"{MULTIQC}/multiqc_report.html")
    conda: "../../envs/preprocessing/multiqc.yaml"
    shell:
        """
        multiqc results logs -o {output}
        """