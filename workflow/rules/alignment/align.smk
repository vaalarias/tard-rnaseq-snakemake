# workflow/rules/alignment/align.smk
if ALIGNER == "star":

    rule star_index:
        input:
            fa=GENOME,
            gtf=GTF
        output:
            index=directory(STAR_INDEX)
        params:
            sjdb_overhang=config["params"].get("star_sjdb_overhang", 100)
        threads: 8
        log:
            f"{LOGS}/align/star_index.log"
        conda:
            "../../envs/alignment/star.yaml"
        shell:
            """
            mkdir -p {output.index} $(dirname {log})

            STAR --runMode genomeGenerate \
                --genomeDir {output.index} \
                --genomeFastaFiles {input.fa:q} \
                --sjdbGTFfile {input.gtf:q} \
                --sjdbOverhang {params.sjdb_overhang} \
                --runThreadN {threads} > {log} 2>&1
            """


    rule star_align:
        input:
            r1=f"{TRIMMED}/{{sample}}_R1_val_1.fq.gz",
            r2=f"{TRIMMED}/{{sample}}_R2_val_2.fq.gz",
            index=STAR_INDEX
        output:
            bam=f"{ALIGN}/star/{{sample}}.bam",
            bai=f"{ALIGN}/star/{{sample}}.bam.bai",
            unmapped_r1=f"{ALIGN}/star/unmapped/{{sample}}_R1_unmapped.fq",
            unmapped_r2=f"{ALIGN}/star/unmapped/{{sample}}_R2_unmapped.fq",
            log_final=f"{ALIGN}/star/logs/{{sample}}_Log.final.out",
            log_out=f"{ALIGN}/star/logs/{{sample}}_Log.out",
            sj=f"{ALIGN}/star/logs/{{sample}}_SJ.out.tab"
        params:
            opts=config["params"]["star"],
            prefix=lambda wildcards: f"{ALIGN}/star/tmp/{wildcards.sample}."
        threads: 8
        log:
            f"{LOGS}/align/{{sample}}_star.log"
        conda:
            "../../envs/alignment/star.yaml"
        shell:
            """
            mkdir -p {ALIGN}/star {ALIGN}/star/unmapped {ALIGN}/star/logs {ALIGN}/star/tmp $(dirname {log})

            STAR \
                --genomeDir {input.index:q} \
                --readFilesIn {input.r1:q} {input.r2:q} \
                --readFilesCommand zcat \
                --runThreadN {threads} \
                {params.opts} \
                --outFileNamePrefix {params.prefix} > {log} 2>&1

            mv {params.prefix}Aligned.sortedByCoord.out.bam {output.bam}
            mv {params.prefix}Unmapped.out.mate1 {output.unmapped_r1}
            mv {params.prefix}Unmapped.out.mate2 {output.unmapped_r2}
            mv {params.prefix}Log.final.out {output.log_final}
            mv {params.prefix}Log.out {output.log_out}
            mv {params.prefix}SJ.out.tab {output.sj}

            samtools index {output.bam} {output.bai}
            """

if ALIGNER == "hisat":

    rule hisat2_index:
        input:
            fa=GENOME
        output:
            index=expand(f"{HISAT_INDEX}/genome.{{i}}.ht2", i=range(1, 9))
        threads: 8
        log:
            f"{LOGS}/align/hisat2_index.log"
        conda:
            "../../envs/alignment/hisat.yaml"
        shell:
            """
            mkdir -p {HISAT_INDEX} $(dirname {log})

            hisat2-build \
                -p {threads} \
                {input.fa:q} \
                {HISAT_INDEX}/genome > {log} 2>&1
            """


    rule hisat2_align:
        input:
            index=rules.hisat2_index.output.index,
            r1=f"{TRIMMED}/{{sample}}_R1_val_1.fq.gz",
            r2=f"{TRIMMED}/{{sample}}_R2_val_2.fq.gz"
        output:
            bam=f"{ALIGN}/hisat/{{sample}}.bam",
            bai=f"{ALIGN}/hisat/{{sample}}.bam.bai",
            unmapped_r1=f"{ALIGN}/hisat/unmapped/{{sample}}_R1_unmapped.fq.gz",
            unmapped_r2=f"{ALIGN}/hisat/unmapped/{{sample}}_R2_unmapped.fq.gz"
        params:
            opts=config["params"].get("hisat2", ""),
            index_prefix=f"{HISAT_INDEX}/genome",
            unmapped_prefix=lambda wildcards: f"{ALIGN}/hisat/unmapped/{wildcards.sample}_unmapped.fq.gz"
        threads: 8
        log:
            f"{LOGS}/align/{{sample}}_hisat2.log"
        conda:
            "../../envs/alignment/hisat.yaml"
        shell:
            """
            mkdir -p {ALIGN}/hisat {ALIGN}/hisat/unmapped $(dirname {log})

            hisat2 \
                -x {params.index_prefix:q} \
                -1 {input.r1:q} \
                -2 {input.r2:q} \
                --threads {threads} \
                {params.opts} \
                --un-conc-gz {params.unmapped_prefix} 2> {log} \
            | samtools sort -@ {threads} -o {output.bam}

            mv {ALIGN}/hisat/unmapped/{wildcards.sample}_unmapped.1.fq.gz {output.unmapped_r1}
            mv {ALIGN}/hisat/unmapped/{wildcards.sample}_unmapped.2.fq.gz {output.unmapped_r2}

            samtools index {output.bam} {output.bai}
            """
            
if SALMON:

    rule salmon_index:
        input:
            gentrome=GENTROME,
            decoy=DECOY
        output:
            index=directory(SALMON_INDEX)
        params:
            k=config["params"].get("salmon_k", 31)
        threads: 8
        log:
            f"{LOGS}/align/salmon_index.log"
        conda:
            "../../envs/alignment/salmon.yaml"
        shell:
            """
            mkdir -p {output.index} $(dirname {log})

            salmon index \
                -t {input.gentrome:q} \
                -d {input.decoy:q} \
                -p {threads} \
                -i {output.index} \
                -k {params.k} > {log} 2>&1
            """


    rule salmon_quant:
        input:
            index=SALMON_INDEX,
            r1=f"{TRIMMED}/{{sample}}_R1_val_1.fq.gz",
            r2=f"{TRIMMED}/{{sample}}_R2_val_2.fq.gz"
        output:
            quant=directory(f"{ALIGN}/salmon/{{sample}}"),
            sf=f"{ALIGN}/salmon/{{sample}}/quant.sf"
        params:
            opts=config["params"].get("salmon", "")
        threads: 8
        log:
            f"{LOGS}/align/{{sample}}_salmon.log"
        conda:
            "../../envs/alignment/salmon.yaml"
        shell:
            """
            mkdir -p {ALIGN}/salmon $(dirname {log})

            salmon quant \
                -i {input.index:q} \
                -l A \
                -1 {input.r1:q} \
                -2 {input.r2:q} \
                -p {threads} \
                {params.opts} \
                -o {output.quant} > {log} 2>&1
            """