import argparse
import io
from pathlib import Path

import pandas as pd


parser = argparse.ArgumentParser(
    usage="%(prog)s INPUT_FILE OUTPUT_DIR",  # Edit this usage line as needed.
    description="Extract tree and cluster sections into TSV files.",
)
parser.add_argument(
    "input_file",
    type=Path,
    help="Input text file.",
)
parser.add_argument(
    "output_dir",
    type=Path,
    help="Existing directory for the output TSV files.",
)
args = parser.parse_args()

if not args.output_dir.is_dir():
    parser.error("OUTPUT_DIR must be an existing directory.")

with args.input_file.open("r", encoding="utf-8") as f:
    sections = f.read().split("\n\n")

outdir = args.output_dir
for section in sections:
    if section.startswith("PARENT_NODE"):
        tree = pd.read_csv(io.StringIO(section), sep="\t")
        tree = tree.rename(columns={"PARENT_NODE": "parent", "CHILD_NODE": "child"})
        tree = tree.sort_values("child")
        tree.to_csv(outdir / "fixed_B.tsv", sep="\t", index=False)
    elif section.startswith("CLUSTER"):
        delta = pd.read_csv(io.StringIO(section), sep="\t")
        delta = delta.rename(columns={"CLUSTER": "cluster", "NODE_OF_CLUSTER": "node"})
        delta.to_csv(outdir / "fixed_delta.tsv", sep="\t", index=False)
