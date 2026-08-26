#!/usr/bin/env python3

from collections import Counter
from pathlib import Path
import gzip
import random


def open_fastq(path):
    path = str(path)

    if path.endswith(".gz"):
        return gzip.open(path, "rt")

    return open(path, "rt", encoding="utf-8")


def read_fastq_sequences(path):
    with open_fastq(path) as handle:
        while True:
            header = handle.readline()

            if not header:
                break

            sequence = handle.readline().strip()
            plus = handle.readline()
            quality = handle.readline()

            if not sequence or not plus or not quality:
                raise ValueError(
                    f"Incomplete FASTQ record in {path}"
                )

            yield sequence.upper()


sample = str(snakemake.wildcards.sample)

r1 = str(snakemake.input["r1"])
r2 = str(snakemake.input["r2"])

for fastq in (r1, r2):
    if not Path(fastq).is_file():
        raise FileNotFoundError(
            f"Unmapped FASTQ does not exist: {fastq}"
        )

    if Path(fastq).stat().st_size == 0:
        raise RuntimeError(
            f"Unmapped FASTQ is empty: {fastq}"
        )
               
summary_output = Path(
    str(snakemake.output["summary"])
)

fasta_output = Path(
    str(snakemake.output["fasta"])
)

log_output = Path(
    str(snakemake.log[0])
)

head_n = int(
    snakemake.params["head_n"]
)

fraction = float(
    snakemake.params["fraction"]
)

min_length = int(
    snakemake.params["min_length"]
)

max_n_fraction = float(
    snakemake.params["max_n_fraction"]
)

seed = int(
    snakemake.params["seed"]
)

if not 0 < fraction <= 1:
    raise ValueError(
        "subsample_fraction must be greater than 0 and at most 1"
    )

summary_output.parent.mkdir(
    parents=True,
    exist_ok=True
)

fasta_output.parent.mkdir(
    parents=True,
    exist_ok=True
)

log_output.parent.mkdir(
    parents=True,
    exist_ok=True
)

rng = random.Random(seed)

sequence_counts = Counter()

total_reads = 0
subsampled_reads = 0
accepted_reads = 0
short_reads = 0
ambiguous_reads = 0

for fastq in (r1, r2):
    for sequence in read_fastq_sequences(fastq):
        total_reads += 1

        if rng.random() > fraction:
            continue

        subsampled_reads += 1

        if len(sequence) < min_length:
            short_reads += 1
            continue

        n_fraction = (
            sequence.count("N") /
            len(sequence)
        )

        if n_fraction > max_n_fraction:
            ambiguous_reads += 1
            continue

        sequence_counts[sequence] += 1
        accepted_reads += 1

top_sequences = sorted(
    sequence_counts.items(),
    key=lambda item: (
        -item[1],
        item[0]
    )
)[:head_n]

with summary_output.open(
    "w",
    encoding="utf-8"
) as handle:
    handle.write(
        "sample\trank\tcount\tlength\tsequence\n"
    )

    for rank, (sequence, count) in enumerate(
        top_sequences,
        start=1
    ):
        handle.write(
            f"{sample}\t"
            f"{rank}\t"
            f"{count}\t"
            f"{len(sequence)}\t"
            f"{sequence}\n"
        )

with fasta_output.open(
    "w",
    encoding="utf-8"
) as handle:
    for rank, (sequence, count) in enumerate(
        top_sequences,
        start=1
    ):
        query_id = (
            f"{sample}"
            f"|rank={rank}"
            f"|count={count}"
            f"|length={len(sequence)}"
        )

        handle.write(
            f">{query_id}\n{sequence}\n"
        )

with log_output.open(
    "w",
    encoding="utf-8"
) as handle:
    handle.write(
        f"sample\t{sample}\n"
        f"r1\t{r1}\n"
        f"r2\t{r2}\n"
        f"subsample_fraction\t{fraction}\n"
        f"total_reads\t{total_reads}\n"
        f"subsampled_reads\t{subsampled_reads}\n"
        f"accepted_reads\t{accepted_reads}\n"
        f"short_reads\t{short_reads}\n"
        f"ambiguous_reads\t{ambiguous_reads}\n"
        f"unique_sequences\t{len(sequence_counts)}\n"
        f"reported_sequences\t{len(top_sequences)}\n"
    )

if not top_sequences:
    raise RuntimeError(
        f"No acceptable unmapped sequences found for {sample}"
    )