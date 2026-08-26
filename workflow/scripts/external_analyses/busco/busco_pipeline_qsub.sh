#!/usr/bin/env bash
#
# Wrapper to run run_busco_pipeline.sh as an SGE job (qsub).
#
# SUBMIT THE JOB WITH:
#   qsub busco_pipeline_qsub.sh
#
# (the queue is already set below with "#$ -q long"; use -q default
#  instead if your cluster/admin recommends it for this type of job;
#  "long" is usually for jobs that take longer, which fits here since
#  it includes genome download + 5 BUSCO runs)
#
# ---------------------------------------------------------------------

#$ -N busco_tardigrade
#$ -cwd
#$ -j y
#$ -o busco_pipeline.log
#$ -pe smp 8
#$ -q long
#$ -S /bin/bash

# -N     : job name (shown in qstat)
# -cwd   : run from the current directory (where you ran qsub), not $HOME
# -j y   : merge stdout and stderr into a single log
# -o     : log file
# -pe smp 8 : request 8 parallel cores (same as THREADS in the script)
# -q long   : queue to use (see "qstat -g c" for available queues)
# -S     : shell to run the script with

set -euo pipefail

echo "=== SGE job started: $(date) ==="
echo "Running on node: $(hostname)"
echo "Working directory: $(pwd)"
echo

# Needed so "conda activate" works inside a non-interactive job
CONDA_BASE="$(conda info --base)"
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# Run the main pipeline
bash run_busco_pipeline.sh

echo
echo "=== SGE job finished: $(date) ==="