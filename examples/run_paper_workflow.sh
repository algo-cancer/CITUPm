#!/usr/bin/env bash
# Run from any directory. PYTHON may select an activated environment's interpreter.
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
python_bin="${PYTHON:-python}"
output_dir="${1:-$repo_dir/results/example}"
mkdir -p "$output_dir/fit" "$output_dir/joint"
common=(--frequenciesFile "$repo_dir/examples/frequencies.tsv"
        --objective cluster --clusterWeight equal --clip-observed-at 1.0
        --presence-threshold 0.025 --max-root-out-degree 1 --numSolutions 1
        --gurobi-threads 1)

# Stage 1: determine E_hat without migration constraints.
"$python_bin" "$repo_dir/citupm.py" "${common[@]}" -o "$output_dir/fit/example"

# Preserve the full-precision reference and use the first tree as a partial MIP start.
reference_score="$("$python_bin" - "$output_dir/fit/example.tree_0.txt" "$output_dir/fit" <<'PY'
from pathlib import Path
import sys
sections = Path(sys.argv[1]).read_text().split('\n\n')
lines = sections[0].splitlines()
header = dict(zip(lines[::2], lines[1::2]))
if header['FIT_OPTIMAL'] != 'YES':
    raise SystemExit('Stage 1 was not proven optimal; inspect the solver log.')
for section, filename in [('PARENT_NODE\tCHILD_NODE', 'B.tsv'),
                          ('CLUSTER\tNODE_OF_CLUSTER', 'delta.tsv')]:
    rows = next(block.splitlines()[1:] for block in sections if block.startswith(section))
    Path(sys.argv[2], filename).write_text('\n'.join(rows) + '\n')
print(header['TREE_SCORE_NO_REWARD'])
PY
)"

# Stages 2 and 3: minimize r under E <= (sqrt(E_hat) + epsilon)^2,
# then fix r and minimize E. Recursive admissibility is always active.
"$python_bin" "$repo_dir/citupm.py" "${common[@]}" -o "$output_dir/joint/example" \
    --sites-file "$repo_dir/examples/sites.tsv" --metastatic-sites metastasis \
    --minimize-migrations-at-score "$reference_score" --migration-score-epsilon 0.05 \
    --warm-start-first-phase-b "$output_dir/fit/B.tsv" \
    --warm-start-first-phase-delta "$output_dir/fit/delta.tsv" \
    --max-refined-solutions 1

"$python_bin" "$repo_dir/plot_with_gene_names.py" "$output_dir/joint" \
    --sample-sites "$repo_dir/examples/sites.tsv" --with-site-transitions --with-node-info
"$python_bin" "$repo_dir/plot_comigration_graphs.py" "$output_dir/joint" \
    --sites-file "$repo_dir/examples/sites.tsv"
