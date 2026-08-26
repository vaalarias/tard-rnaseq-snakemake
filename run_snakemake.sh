#!/bin/bash
#
# Run the RNA-seq Snakemake pipeline on the cluster (SGE/SLURM).
# use profiles/sge or profiles/slurm respectively
# Run this script from the repository root.

NOTIFY_EMAIL="Write here some email :)"

notify_when_finished() {
    exit_status=$?

    trap - EXIT

    if [ "$exit_status" -eq 0 ]; then
        subject="Snakemake completed successfully"
        result="SUCCESS"
    else
        subject="Snakemake failed"
        result="FAILED"
    fi

    {
        printf 'Workflow status: %s\n' "$result"
        printf 'Exit status: %s\n' "$exit_status"
        printf 'Host: %s\n' "$(hostname)"
        printf 'Directory: %s\n' "$(pwd)"
        printf 'Finished: %s\n' "$(date)"
        printf 'Targets/arguments: %s\n' "$*"
    } |
        mail -s "$subject" "$NOTIFY_EMAIL" ||
        true

    exit "$exit_status"
}

trap 'notify_when_finished "$@"' EXIT

set -euo pipefail

WORK_CACHE="/export/space3/users/vjarias/.snakemake_cache"

export TMPDIR="${WORK_CACHE}/tmp"
export CONDA_PKGS_DIRS="${WORK_CACHE}/conda_pkgs"
export XDG_CACHE_HOME="${WORK_CACHE}/xdg"

mkdir -p \
    "${TMPDIR}" \
    "${CONDA_PKGS_DIRS}" \
    "${XDG_CACHE_HOME}"

echo "TMPDIR: ${TMPDIR}"
echo "Conda package cache: ${CONDA_PKGS_DIRS}"

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

mkdir -p "${TMPDIR}"

echo "Temporary directory set to: ${TMPDIR}"

snakemake \
    "$@" \
    --profile profiles/sge \
    --jobs 10 \
    --use-conda \
    --conda-frontend conda \
    --rerun-incomplete