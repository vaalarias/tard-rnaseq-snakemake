#!/usr/bin/env python3
"""
Snakemake SGE (qsub) submit script.

Reads the job properties Snakemake generates for each job (rule name,
threads, resources) and builds/submits a matching `qsub` command.
Prints ONLY the numeric job ID to stdout on success -- Snakemake
requires that to track the job.
"""

import os
import sys
import subprocess
import re

from snakemake.utils import read_job_properties

jobscript = sys.argv[-1]
job_properties = read_job_properties(jobscript)

rule_name = job_properties.get("rule", "job")
threads = job_properties.get("threads", 1)
resources = job_properties.get("resources", {})

mem_mb = int(resources.get("mem_mb", 4000))
time_min = int(resources.get("time_min", 240))
queue = resources.get("queue", "long")

threads = int(threads)

mem_per_slot_mb = max(1, mem_mb // max(1, threads))

hours, minutes = divmod(time_min, 60)
time_sge = f"{hours:02d}:{minutes:02d}:00"

os.makedirs("logs/sge", exist_ok=True)

omit_h_vmem = int(
    resources.get("omit_h_vmem", 0)
) == 1

log_out = f"logs/sge/{rule_name}.$JOB_ID.out"
log_err = f"logs/sge/{rule_name}.$JOB_ID.err"

qsub_cmd = [
    "qsub",
    "-V",
    "-cwd",
    "-j", "n",
    "-o", log_out,
    "-e", log_err,
    "-N", f"smk_{rule_name}",
    "-pe", "smp", str(threads),
]

if not omit_h_vmem:
    qsub_cmd.extend([
        "-l",
        f"h_vmem={mem_per_slot_mb}M",
    ])

qsub_cmd.extend([
    "-l", f"h_rt={time_sge}",
    "-q", queue,
    "-terse",
    jobscript,
])


result = subprocess.run(
    qsub_cmd,
    capture_output=True,
    text=True
)

if result.returncode != 0:
    sys.stderr.write(
        "SGE submission failed.\n"
        f"Command: {' '.join(qsub_cmd)}\n"
        f"Return code: {result.returncode}\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}\n"
    )
    sys.exit(result.returncode)
    
# With -terse, qsub prints just the job ID (possibly with a trailing
# newline / task range for array jobs -- not used here, so this is safe)
job_id_match = re.search(r"\d+", result.stdout)
if not job_id_match:
    sys.stderr.write(f"Could not parse job ID from qsub output: {result.stdout}\n")
    sys.exit(1)

print(job_id_match.group())
