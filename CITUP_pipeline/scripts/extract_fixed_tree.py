import sys
import pandas as pd
import io
from pathlib import Path

with open(sys.argv[1], "r", encoding="utf-8") as f:
    sections = f.read().split("\n\n")

outdir = Path(sys.argv[2])
for section in sections:
    if section.startswith('PARENT_NODE'):
        tree = pd.read_csv(io.StringIO(section), sep="\t")
        tree = tree.rename(columns={"PARENT_NODE": 'parent', "CHILD_NODE": 'child'})
        tree = tree.sort_values('child')
        tree.to_csv(outdir / 'fixed_B.tsv', sep="\t", index=False)
    elif section.startswith("CLUSTER"):
        delta = pd.read_csv(io.StringIO(section), sep="\t")
        delta = delta.rename(columns={"CLUSTER": 'cluster', "NODE_OF_CLUSTER": 'node'})
        delta.to_csv(outdir / 'fixed_delta.tsv', sep="\t", index=False)
