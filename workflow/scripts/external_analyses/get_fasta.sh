#!/bin/bash
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --mem=50G
#SBATCH --cpus-per-task=8
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1

#SBATCH --job-name=getfasta
#SBATCH --output=%x_%j_output.txt
#SBATCH --error=%x_%j_error.txt


###stringtie /lustre/groups/ife/workspace/schneider/valentina.ojeda/rnaseq/explore/dataset2_results2/alignments/star/sorted_ECO-3.bam -G /lustre/groups/ife/workspace/schneider/valentina.ojeda/rnaseq/snakemake1/workflow/scripts/GCF_019649055.1_Prichtersi_v1.0_genomic.gff -o ECO-3.gtf -p 8 -e 
###stringtie /lustre/groups/ife/workspace/schneider/valentina.ojeda/rnaseq/explore/dataset2_results2/alignments/star/sorted_GCB-4.bam -G /lustre/groups/ife/workspace/schneider/valentina.ojeda/rnaseq/snakemake1/workflow/scripts/GCF_019649055.1_Prichtersi_v1.0_genomic.gff -o GCB-4.gtf -p 8 -e 
bedtools getfasta -fi /lustre/groups/ife/workspace/schneider/valentina.ojeda/rnaseq/snakemake1/config/ref/GCF_019649055.1_Prichtersi_v1.0_genomic.fna -bed ECO-3.gtf -fo ECO-3_transcripts.fasta -name
bedtools getfasta -fi /lustre/groups/ife/workspace/schneider/valentina.ojeda/rnaseq/snakemake1/config/ref/GCF_019649055.1_Prichtersi_v1.0_genomic.fna -bed GCB-4.gtf -fo GCB-4_transcripts.fasta -name