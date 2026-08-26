#!/usr/bin/env python3

from pathlib import Path
import hashlib
import os
import shutil
import subprocess
import time


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

chunk_size = int(
    snakemake.params["chunk_size"]
)

max_retries = int(
    snakemake.params["max_retries"]
)

retry_wait_seconds = int(
    snakemake.params["retry_wait_seconds"]
)

attempt_timeout_seconds = int(
    snakemake.params["attempt_timeout_seconds"]
)

output_path.parent.mkdir(
    parents=True,
    exist_ok=True
)

log_path.parent.mkdir(
    parents=True,
    exist_ok=True
)

chunk_directory = (
    output_path.parent /
    "blast_chunks"
)

chunk_directory.mkdir(
    parents=True,
    exist_ok=True
)

if not query_path.is_file():
    raise FileNotFoundError(
        f"BLAST query FASTA does not exist: {query_path}"
    )

if query_path.stat().st_size == 0:
    raise RuntimeError(
        f"BLAST query FASTA is empty: {query_path}"
    )

if chunk_size < 1:
    raise ValueError(
        "BLAST chunk size must be at least 1"
    )

if max_retries < 1:
    raise ValueError(
        "BLAST max_retries must be at least 1"
    )

if attempt_timeout_seconds < 1:
    raise ValueError(
        "BLAST attempt timeout must be at least 1 second"
    )


def read_fasta_records(path):
    header = None
    sequence_lines = []

    with path.open(
        "r",
        encoding="utf-8"
    ) as handle:
        for line in handle:
            line = line.rstrip("\n")

            if line.startswith(">"):
                if header is not None:
                    yield header, sequence_lines

                header = line
                sequence_lines = []

            else:
                if header is None:
                    raise ValueError(
                        f"Sequence found before FASTA "
                        f"header in {path}"
                    )

                if line:
                    sequence_lines.append(line)

    if header is not None:
        yield header, sequence_lines


def format_fasta_records(records):
    lines = []

    for header, sequence_lines in records:
        lines.append(header)
        lines.extend(sequence_lines)

    return "\n".join(lines) + "\n"


def chunk_identifier(fasta_text):
    parameter_text = (
        f"database={database}\n"
        f"task={task}\n"
        f"evalue={evalue}\n"
        f"max_target_seqs={max_target_seqs}\n"
    )

    digest = hashlib.sha256(
        (
            parameter_text +
            fasta_text
        ).encode("utf-8")
    ).hexdigest()

    return digest[:12]


records = list(
    read_fasta_records(query_path)
)

if not records:
    raise RuntimeError(
        f"No FASTA records found in {query_path}"
    )

outfmt = (
    "6 qseqid saccver pident length mismatch gapopen "
    "qstart qend sstart send evalue bitscore "
    "staxids sscinames sskingdoms stitle"
)

chunk_jobs = []

for chunk_start in range(
    0,
    len(records),
    chunk_size
):
    chunk_number = (
        chunk_start // chunk_size
    ) + 1

    chunk_records = records[
        chunk_start:
        chunk_start + chunk_size
    ]

    fasta_text = format_fasta_records(
        chunk_records
    )

    identifier = chunk_identifier(
        fasta_text
    )

    prefix = (
        chunk_directory /
        f"chunk_{chunk_number:03d}_{identifier}"
    )

    chunk_path = Path(
        str(prefix) + ".fa"
    )

    chunk_output = Path(
        str(prefix) + ".blast.tsv"
    )

    completion_marker = Path(
        str(prefix) + ".done"
    )

    chunk_path.write_text(
        fasta_text,
        encoding="utf-8"
    )

    chunk_jobs.append({
        "number": chunk_number,
        "query": chunk_path,
        "output": chunk_output,
        "done": completion_marker,
    })

temporary_final = output_path.with_suffix(
    output_path.suffix + ".tmp"
)

try:
    with log_path.open(
        "w",
        encoding="utf-8"
    ) as log_handle:
        log_handle.write(
            f"query\t{query_path}\n"
            f"sequences\t{len(records)}\n"
            f"chunk_size\t{chunk_size}\n"
            f"chunks\t{len(chunk_jobs)}\n"
            f"max_retries\t{max_retries}\n"
            f"attempt_timeout_seconds\t"
            f"{attempt_timeout_seconds}\n"
        )

        log_handle.flush()

        with temporary_final.open(
            "w",
            encoding="utf-8"
        ) as combined_handle:
            for chunk_job in chunk_jobs:
                chunk_index = chunk_job["number"]
                chunk_path = chunk_job["query"]
                chunk_output = chunk_job["output"]
                completion_marker = chunk_job["done"]

                if (
                    completion_marker.exists() and
                    chunk_output.exists()
                ):
                    log_handle.write(
                        "\n"
                        f"Reusing completed chunk "
                        f"{chunk_index}/{len(chunk_jobs)}\n"
                    )

                    log_handle.flush()

                else:
                    success = False

                    for attempt in range(
                        1,
                        max_retries + 1
                    ):
                        temporary_chunk_output = Path(
                            str(chunk_output) + ".tmp"
                        )

                        if temporary_chunk_output.exists():
                            temporary_chunk_output.unlink()

                        log_handle.write(
                            "\n"
                            f"Running chunk "
                            f"{chunk_index}/{len(chunk_jobs)}, "
                            f"attempt "
                            f"{attempt}/{max_retries}\n"
                        )

                        log_handle.flush()

                        command = [
                            "blastn",
                            "-query",
                            str(chunk_path),
                            "-db",
                            database,
                            "-task",
                            task,
                            "-remote",
                            "-evalue",
                            evalue,
                            "-max_target_seqs",
                            max_target_seqs,
                            "-outfmt",
                            outfmt,
                            "-out",
                            str(temporary_chunk_output),
                        ]

                        try:
                            completed = subprocess.run(
                                command,
                                stdout=log_handle,
                                stderr=subprocess.STDOUT,
                                check=False,
                                timeout=(
                                    attempt_timeout_seconds
                                ),
                            )

                            return_code = (
                                completed.returncode
                            )

                        except subprocess.TimeoutExpired:
                            return_code = None

                            log_handle.write(
                                f"Chunk {chunk_index}, "
                                f"attempt {attempt} exceeded "
                                f"{attempt_timeout_seconds} "
                                "seconds\n"
                            )

                            log_handle.flush()

                        if return_code == 0:
                            if not temporary_chunk_output.exists():
                                temporary_chunk_output.touch()

                            os.replace(
                                temporary_chunk_output,
                                chunk_output
                            )

                            completion_marker.write_text(
                                "completed\n",
                                encoding="utf-8"
                            )

                            success = True

                            log_handle.write(
                                f"Chunk {chunk_index} "
                                "completed\n"
                            )

                            log_handle.flush()
                            break

                        if temporary_chunk_output.exists():
                            temporary_chunk_output.unlink()

                        log_handle.write(
                            f"Chunk {chunk_index} failed"
                        )

                        if return_code is not None:
                            log_handle.write(
                                f" with exit code "
                                f"{return_code}"
                            )

                        log_handle.write("\n")
                        log_handle.flush()

                        if attempt < max_retries:
                            wait_seconds = (
                                retry_wait_seconds *
                                attempt
                            )

                            log_handle.write(
                                f"Waiting {wait_seconds} "
                                "seconds before retry\n"
                            )

                            log_handle.flush()

                            time.sleep(
                                wait_seconds
                            )

                    if not success:
                        raise RuntimeError(
                            f"Remote BLAST failed after "
                            f"{max_retries} attempts for "
                            f"{chunk_path.name}"
                        )

                with chunk_output.open(
                    "r",
                    encoding="utf-8"
                ) as chunk_handle:
                    shutil.copyfileobj(
                        chunk_handle,
                        combined_handle
                    )

        os.replace(
            temporary_final,
            output_path
        )

        log_handle.write(
            "\nRemote BLAST completed successfully\n"
        )

        log_handle.flush()

finally:
    if temporary_final.exists():
        temporary_final.unlink()