rule featurecounts:
    input:
        bam=f"{ALIGN}/star/{{sample}}.bam"
    output:
        counts=f"{COUNTS}/{{sample}}_counts.txt"
    threads: 4
    params:
        gtf=GTF,
        opts=config["params"]["featurecounts"]
    log:
        f"{LOGS}/featurecounts/{{sample}}.log"
    conda:
        "../../envs/alignment/featurecounts.yaml"
    shell:
        """
        mkdir -p {COUNTS} $(dirname {log})

        featureCounts \
            -T {threads} \
            -a {params.gtf:q} \
            {params.opts} \
            -o {output.counts:q} \
            {input.bam:q} > {log} 2>&1
        """