#!/usr/bin/env python3

from pathlib import Path
import os
import subprocess


query_path = Path(
    str(snakemake.input["fasta"])
)

output_path = Path(
    str(snakemake.output["blast"])
)

log_path = Path(
    str(snakemake.log[0])
)

database = str(
    snakemake.params["database"]
)

task = str(
    snakemake.params["task"]
)

evalue = str(
    snakemake.params["evalue"]
)

max_target_seqs = str(
    snakemake.params["max_target_seqs"]
)

threads = int(
    snakemake.threads
)

if not query_path.is_file():
    raise FileNotFoundError(
        f"BLAST query FASTA does not exist: {query_path}"
    )

if query_path.stat().st_size == 0:
    raise RuntimeError(
        f"BLAST query FASTA is empty: {query_path}"
    )

output_path.parent.mkdir(
    parents=True,
    exist_ok=True
)

log_path.parent.mkdir(
    parents=True,
    exist_ok=True
)

temporary_output = output_path.with_suffix(
    output_path.suffix + ".tmp"
)

outfmt = (
    "6 qseqid saccver pident length mismatch gapopen "
    "qstart qend sstart send evalue bitscore "
    "staxids sscinames sskingdoms stitle"
)

command = [
    "blastn",
    "-query",
    str(query_path),
    "-db",
    database,
    "-task",
    task,
    "-evalue",
    evalue,
    "-max_target_seqs",
    max_target_seqs,
    "-num_threads",
    str(threads),
    "-outfmt",
    outfmt,
    "-out",
    str(temporary_output),
]

try:
    if temporary_output.exists():
        temporary_output.unlink()

    with log_path.open(
        "w",
        encoding="utf-8"
    ) as log_handle:
        log_handle.write(
            f"query\t{query_path}\n"
            f"database\t{database}\n"
            f"task\t{task}\n"
            f"evalue\t{evalue}\n"
            f"max_target_seqs\t{max_target_seqs}\n"
            f"threads\t{threads}\n"
            "mode\tlocal\n\n"
        )

        log_handle.flush()

        completed = subprocess.run(
            command,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            check=False,
        )

        if completed.returncode != 0:
            raise RuntimeError(
                "Local BLAST failed with exit code "
                f"{completed.returncode}. See {log_path}"
            )

        if not temporary_output.exists():
            raise RuntimeError(
                "BLAST finished without creating its output."
            )

        os.replace(
            temporary_output,
            output_path
        )

        log_handle.write(
            "\nLocal BLAST completed successfully\n"
        )

finally:
    if temporary_output.exists():
        temporary_output.unlink()