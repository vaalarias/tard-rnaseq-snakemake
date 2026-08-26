#!/usr/bin/env python3
"""
Snakemake SGE job status checker.

Snakemake calls this with the job ID returned by sge-submit.py, and
expects exactly one of these printed to stdout: "running", "success",
or "failed".
"""

import sys
import subprocess

jobid = sys.argv[1]

# 1) Still in the queue or running?
qstat = subprocess.run(
    ["qstat", "-j", jobid], capture_output=True, text=True
)
if qstat.returncode == 0:
    # qstat -j succeeds while the job is queued (qw) or running (r)
    print("running")
    sys.exit(0)

# 2) Job is no longer in the queue -- check its exit status via qacct.
#    qacct can lag a few seconds behind job completion, so retry a
#    couple of times before giving up.
import time

for _ in range(5):
    qacct = subprocess.run(
        ["qacct", "-j", jobid], capture_output=True, text=True
    )
    if qacct.returncode == 0:
        exit_status = None
        for line in qacct.stdout.splitlines():
            if line.startswith("exit_status"):
                exit_status = line.split()[-1].strip()
                break
        if exit_status == "0":
            print("success")
        else:
            print("failed")
        sys.exit(0)
    time.sleep(3)

# If we still can't find it in qacct, report failed rather than hanging
# Snakemake indefinitely.
print("failed")
