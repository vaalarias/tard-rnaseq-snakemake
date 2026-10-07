# workflow/scripts/preprocessing/summarize_trial_mapping.py

from pathlib import Path


TRIAL_DIR = Path(snakemake.params.trial_dir)
GENOMES = list(snakemake.params.genomes)
SAMPLES = list(snakemake.params.samples)

OUTPUT = Path(snakemake.output.summary)


# STAR metrics to extract from Log.final.out
METRICS = {
    "Number of input reads": "input_reads",

    "Uniquely mapped reads number": "uniquely_mapped_reads",
    "Uniquely mapped reads %": "uniquely_mapped_pct",

    "Number of reads mapped to multiple loci": "multimapped_reads",
    "% of reads mapped to multiple loci": "multimapped_pct",

    "Number of reads mapped to too many loci": "too_many_loci_reads",
    "% of reads mapped to too many loci": "too_many_loci_pct",

    "Number of reads unmapped: too many mismatches":
        "unmapped_mismatches_reads",
    "% of reads unmapped: too many mismatches":
        "unmapped_mismatches_pct",

    "Number of reads unmapped: too short":
        "unmapped_too_short_reads",
    "% of reads unmapped: too short":
        "unmapped_too_short_pct",

    "Number of reads unmapped: other":
        "unmapped_other_reads",
    "% of reads unmapped: other":
        "unmapped_other_pct",
}


COLUMNS = [
    "genome",
    "sample",
    "input_reads",
    "uniquely_mapped_reads",
    "uniquely_mapped_pct",
    "multimapped_reads",
    "multimapped_pct",
    "too_many_loci_reads",
    "too_many_loci_pct",
    "unmapped_mismatches_reads",
    "unmapped_mismatches_pct",
    "unmapped_too_short_reads",
    "unmapped_too_short_pct",
    "unmapped_other_reads",
    "unmapped_other_pct",
]


def parse_star_log(log_file):

    values = {}

    with open(log_file) as handle:

        for line in handle:

            if "|" not in line:
                continue

            key, value = line.split("|", 1)

            key = key.strip()
            value = value.strip()

            if key in METRICS:
                values[METRICS[key]] = value.rstrip("%")

    return values


OUTPUT.parent.mkdir(parents=True, exist_ok=True)


with open(OUTPUT, "w") as out:

    out.write("\t".join(COLUMNS) + "\n")

    for genome in GENOMES:

        for sample in SAMPLES:

            log_file = (
                TRIAL_DIR
                / genome
                / f"{sample}_Log.final.out"
            )

            values = parse_star_log(log_file)

            row = [
                genome,
                sample,
                *[
                    values.get(column, "NA")
                    for column in COLUMNS[2:]
                ],
            ]

            out.write("\t".join(row) + "\n")