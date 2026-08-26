#!/usr/bin/env python3

from pathlib import Path
import re

import pandas as pd


blast_path = Path(
    str(snakemake.input["blast"])
)

fasta_path = Path(
    str(snakemake.input["fasta"])
)

best_hits_output = Path(
    str(snakemake.output["best_hits"])
)

taxon_summary_output = Path(
    str(snakemake.output["taxon_summary"])
)

best_hits_output.parent.mkdir(
    parents=True,
    exist_ok=True
)

taxon_summary_output.parent.mkdir(
    parents=True,
    exist_ok=True
)

blast_columns = [
    "query_id",
    "subject_accession",
    "percent_identity",
    "alignment_length",
    "mismatches",
    "gap_opens",
    "query_start",
    "query_end",
    "subject_start",
    "subject_end",
    "evalue",
    "bit_score",
    "taxid",
    "scientific_name",
    "kingdom",
    "subject_title",
]


def parse_query_id(query_id):
    sample_match = re.match(
        r"^([^|]+)",
        query_id
    )

    rank_match = re.search(
        r"\|rank=(\d+)",
        query_id
    )

    count_match = re.search(
        r"\|count=(\d+)",
        query_id
    )

    length_match = re.search(
        r"\|length=(\d+)",
        query_id
    )

    return {
        "query_id": query_id,
        "sample": (
            sample_match.group(1)
            if sample_match else ""
        ),
        "rank": (
            int(rank_match.group(1))
            if rank_match else pd.NA
        ),
        "read_count": (
            int(count_match.group(1))
            if count_match else 1
        ),
        "query_length": (
            int(length_match.group(1))
            if length_match else pd.NA
        ),
    }


query_records = []

with fasta_path.open(
    "r",
    encoding="utf-8"
) as handle:
    for line in handle:
        if line.startswith(">"):
            query_id = line[1:].strip().split()[0]

            query_records.append(
                parse_query_id(query_id)
            )

queries = pd.DataFrame(query_records)

if blast_path.stat().st_size == 0:
    blast = pd.DataFrame(
        columns=blast_columns
    )
else:
    blast = pd.read_csv(
        blast_path,
        sep="\t",
        names=blast_columns,
        header=None,
        dtype={
            "query_id": str,
            "subject_accession": str,
            "taxid": str,
            "scientific_name": str,
            "kingdom": str,
            "subject_title": str,
        },
    )

if blast.empty:
    best_hits = queries.copy()

    best_hits["hit_status"] = "no_hit"
    best_hits["subject_accession"] = pd.NA
    best_hits["percent_identity"] = pd.NA
    best_hits["alignment_length"] = pd.NA
    best_hits["evalue"] = pd.NA
    best_hits["bit_score"] = pd.NA
    best_hits["taxid"] = pd.NA
    best_hits["scientific_name"] = "No BLAST hit"
    best_hits["kingdom"] = "No BLAST hit"
    best_hits["subject_title"] = pd.NA
else:
    numeric_columns = [
        "percent_identity",
        "alignment_length",
        "evalue",
        "bit_score",
    ]

    for column in numeric_columns:
        blast[column] = pd.to_numeric(
            blast[column],
            errors="coerce"
        )

    blast = blast.sort_values(
        by=[
            "query_id",
            "evalue",
            "bit_score",
            "percent_identity",
        ],
        ascending=[
            True,
            True,
            False,
            False,
        ],
    )

    best_blast = blast.drop_duplicates(
        subset="query_id",
        keep="first"
    )

    best_hits = queries.merge(
        best_blast,
        on="query_id",
        how="left"
    )

    best_hits["hit_status"] = (
        best_hits["subject_accession"]
        .notna()
        .map({
            True: "hit",
            False: "no_hit"
        })
    )

    best_hits["scientific_name"] = (
        best_hits["scientific_name"]
        .fillna("No BLAST hit")
    )

    best_hits["kingdom"] = (
        best_hits["kingdom"]
        .fillna("No BLAST hit")
    )

best_hits = best_hits.sort_values(
    by=[
        "sample",
        "rank"
    ]
)

best_hits.to_csv(
    best_hits_output,
    sep="\t",
    index=False
)

taxon_summary = (
    best_hits
    .groupby(
        [
            "sample",
            "kingdom",
            "scientific_name",
            "taxid",
        ],
        dropna=False
    )
    .agg(
        unique_queries=(
            "query_id",
            "nunique"
        ),
        weighted_read_count=(
            "read_count",
            "sum"
        ),
        mean_percent_identity=(
            "percent_identity",
            "mean"
        ),
        mean_alignment_length=(
            "alignment_length",
            "mean"
        ),
        best_evalue=(
            "evalue",
            "min"
        ),
    )
    .reset_index()
    .sort_values(
        by=[
            "sample",
            "weighted_read_count"
        ],
        ascending=[
            True,
            False
        ]
    )
)

taxon_summary.to_csv(
    taxon_summary_output,
    sep="\t",
    index=False
)