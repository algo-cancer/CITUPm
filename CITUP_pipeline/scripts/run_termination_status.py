import argparse
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd


parser = argparse.ArgumentParser(
    usage="%(prog)s SCAN_DIR",  # Edit this usage line as needed.
    description="Summarize citup logs and Slurm job status as TSV on stdout.",
)
parser.add_argument(
    "scan_dir",
    type=Path,
    help="Directory to search recursively for citup.log files.",
)
args = parser.parse_args()

if not args.scan_dir.is_dir():
    parser.error("SCAN_DIR must be an existing directory.")

infeasible_regex = re.compile("odel is infeasible")

scandir = args.scan_dir

status = []
for filename in scandir.rglob("citup.log"):
    record = {}
    record["person"] = filename.parent.parent.name
    record["infeasible"] = False
    with open(filename) as fh:
        for i, line in enumerate(fh):
            if i == 0:
                line = line.rstrip()
                record["jobid"] = line[8:]
            if re.search(infeasible_regex, line):
                record["infeasible"] = True
    status.append(record)
status = pd.DataFrame(status)

# Unique job IDs, excluding missing values
jobids = status["jobid"].dropna().astype(str).unique()

# One sacct call for all jobs
result = subprocess.run(
    [
        "sacct",
        "-j",
        ",".join(jobids),
        "--noheader",
        "--parsable2",
        "--allocations",
        "--format=JobIDRaw,State,Elapsed",
    ],
    capture_output=True,
    text=True,
    check=True,
)

# Parse sacct output
records = []
for line in result.stdout.splitlines():
    if line.strip():
        jobid, state, runtime = line.split("|", 2)
        records.append((jobid, state, runtime))

term_status = pd.DataFrame(records, columns=["jobid", "status", "runtime"])
status = status.merge(term_status, how="left", on="jobid")

# Must have failed to be declared infeasible
status.loc[
    status.infeasible & (status.status == "FAILED"), "status"
] = "Infeasible"

status[["person", "jobid", "status", "runtime"]].to_csv(
    sys.stdout, sep="\t", index=False
)
