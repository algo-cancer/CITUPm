# CITUPm

CITUPm jointly infers a clonal tree, clone prevalences, and a migration history from mutation-cluster descendant cell fractions (DCFs) measured across tumor samples. It uses Gurobi to optimize the topology directly, then selects migration refinements without reseeding and with the minimum number of comigrations when such refinements exist.

This repository contains the current inference entry point, `citupm.py`, and the two main plotting scripts. 

## Setup

```bash
git clone https://github.com/algo-cancer/CITUPm.git
cd CITUPm
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Install the Graphviz application separately so `dot` is on `PATH`; the Python `graphviz` package alone does not supply it. For example, use `sudo apt-get install graphviz` on Debian/Ubuntu, `brew install graphviz` on macOS, or your cluster's Graphviz module. See [Graphviz installation](https://graphviz.org/download/).

A Gurobi license is required. Configure your institution's license or set `GRB_LICENSE_FILE` to your license file. See [Gurobi Python installation and licensing](https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python).

```bash
python -c 'import gurobipy, pandas, graphviz; print(gurobipy.gurobi.version())'
dot -V
python citupm.py --help
bash examples/run_paper_workflow.sh
```

The example writes solver outputs and PDFs under `results/example/`. Set `PYTHON=/path/to/python` to select another interpreter, or pass an output directory as the script's first argument.

## Inputs

### Mutation frequencies (`-f`, required)

A whitespace-delimited table, preferably TSV, with a header. The first column contains unique mutation IDs, the last column contains cluster IDs, and every intervening column is a sample's DCF. Do not include the healthy root as a mutation cluster.

```text
mutID   primary_1   primary_2   metastasis_1   clusterID
m1      0.80        0.75        0.70           C1
m2      0.82        0.77        0.72           C1
m3      0.40        0.35        0.30           C2
m4      0.42        0.37        0.32           C2
```

DCFs are proportions, not percentages. Use purity-adjusted or unadjusted values consistently across the whole matrix. CITUPm does not infer purity or perform copy-number correction. `--clip-observed-at 1.0` caps values greater than one before fitting, as in the paper's TRACERx protocol.

For `--objective cluster`, the mean of the mutations in each cluster supplies the observed DCF `f[i,s]`. A matrix with one representative row per cluster is also accepted. With `--clusterWeight equal`, every cluster/sample residual has equal weight; with the script default `size`, a cluster's weight is its mutation count. `sqrt` uses the square root of that count. Keep inputs, clipping, filtering, and weights identical across stages.

### Sample locations (`--sites-file`)

A two-column sample-to-site mapping, with an optional `sample site` header. Every frequency-table sample must appear. Multiple samples may map to one site; presence at a site is supported if any sample from that site meets the presence criterion. Sample and site labels must not contain whitespace.

```text
sample          site
primary_1       primary
primary_2       primary
metastasis_1    metastasis
```

Supply this file for migration inference. Omit it for the unconstrained reference fit, and do not set `--metastatic-sites` or stacking rewards in that reference fit: those options also activate migration inference, assuming one sample per site if no mapping is supplied.

`--metastatic-sites metastasis` optionally prevents that site from being the origin of a root child. Accepts comma-separated and/or space-separated site names. Names must occur in the mapping, and at least one site must remain eligible as an origin.

### Fixed trees and warm starts

Node `0` is the healthy root. With the default fixed cluster mapping, non-root node IDs `1..n` follow cluster first-appearance order after filtering.

- `--fixed-B-file`: complete `parent child` edge list using node IDs; optional header `parent child`.
- `--fixed-delta-file`: complete `cluster node` table for retained clusters; optional header `cluster node`.
- `--warm-start-first-phase-b` and `--warm-start-first-phase-delta`: the same formats, used as partial Gurobi starts for migration minimization. They do not constrain the final tree and require `--minimize-migrations-at-score`.

## Paper workflow

The terminology and equations here follow the supplied RECOMB 2027 manuscript. Its default procedure has three optimization stages:

1. Minimize DCF fit error without migration constraints to obtain `E_hat`.
2. Minimize migration count `r` subject to `E <= (sqrt(E_hat) + epsilon)^2` (main Eq. 9).
3. Fix the selected migration count `r_star` and minimize DCF fit error again. Stages 2 and 3 run in a single invocation by default.

After fitting, the program fixes the inferred topology and prevalences and seeks an acyclic site migration graph with the minimum number of comigrations. `--max-refined-solutions` controls how many such refinements are emitted.

The [complete runnable example](examples/run_paper_workflow.sh) also extracts the stage-1 topology and cluster assignment for warm starts. The commands below show the same procedure for your data. Run from the repository root:

```bash
mkdir -p results/patient/fit results/patient/joint

# Stage 1: no sample-site mapping, hence no migration constraints.
python citupm.py -f /path/to/frequencies.tsv -o results/patient/fit/patient \
  --objective cluster --clusterWeight equal --clip-observed-at 1.0 \
  --presence-threshold 0.025 --max-root-out-degree 1 --numSolutions 1 \
  --fit-time-limit 86400 --gurobi-threads 16

# Read the full-precision DCF fit error; first check that the fit was proven optimal.
reference_score="$(python - <<'PY'
from pathlib import Path
lines = Path('results/patient/fit/patient.tree_0.txt').read_text().split('\n\n')[0].splitlines()
header = dict(zip(lines[::2], lines[1::2]))
if header['FIT_OPTIMAL'] != 'YES':
    raise SystemExit('The reference fit was not proven optimal; inspect its solver log.')
print(header['TREE_SCORE_NO_REWARD'])
PY
)"

# Stages 2 and 3: minimum migrations within the tolerance, then best fit at that count.
python citupm.py -f /path/to/frequencies.tsv -o results/patient/joint/patient \
  --objective cluster --clusterWeight equal --clip-observed-at 1.0 \
  --presence-threshold 0.025 --max-root-out-degree 1 --numSolutions 1 \
  --sites-file /path/to/sites.tsv \
  --minimize-migrations-at-score "$reference_score" --migration-score-epsilon 0.05 \
  --max-refined-solutions 1 --fit-time-limit 86400 --gurobi-threads 16
```

The TRACERx protocol uses equal weights, clipping at 1.0, presence threshold `0.025`, root out-degree `1`, one clonal tree, one refined tree, 16 threads, and a 24-hour fit budget. The simulations use the same approach with a 16-hour fit budget. Supply the primary/metastatic annotations with `--metastatic-sites` when available. To use the paper's warm starts, add the two `--warm-start-first-phase-*` files, as demonstrated in the example script.

All optional reward coefficients and `--migration-weight` are zero by default; leave them zero for the paper's three-stage procedure. The script's generic defaults retain `--clusterWeight size`, `--presence-threshold 0.05`, and up to 20 clonal trees, so supply the paper settings explicitly.

A time-limited incumbent is not a proof of optimality. Check `FIT_OPTIMAL`, `MIGRATION_MINIMUM_PROVEN`, `MIGRATION_SCORE_TIEBREAK_PROVEN`, and refinement status. Stages 2 and 3 share the `--fit-time-limit` budget. Refinement has its own solves; `--refinement-time-limit` sets a limit per refinement solve. `--gurobi-max-time` supplies a default limit for each model, while `--enum-time-limit` caps additional topology enumeration.

### Other inference modes

Fit DCFs at an **exact** migration count, for example to explore the paper's fit/migration tradeoff:

```bash
python citupm.py -f /path/to/frequencies.tsv -o results/fixed_r/patient \
  --sites-file /path/to/sites.tsv --num-migrations 5 \
  --objective cluster --clusterWeight equal --clip-observed-at 1.0 \
  --presence-threshold 0.025 --numSolutions 1
```

Create `results/fixed_r/` first. `--max-migrations` is accepted as the manuscript's spelling of this option; despite that spelling, the inherited formulation fixes the count by equality. `--num-migrations` makes that behavior explicit. It cannot be combined with `--minimize-migrations-at-score`.

`--numSolutions K` enumerates up to `K` distinct topologies ranked by DCF fit error, including ties. It may return fewer when there are no more feasible topologies or a resource limit intervenes. Alternatively, use `--enum-max-trees K` with `--enum-fit-abs` and/or `--enum-fit-rel` for a near-optimal window. `--numSolutions` overrides these window settings. In score-bounded mode, enumeration always preserves the error cap and selected migration count.

## Model notation and admissibility

| Paper symbol | Implementation | Meaning |
| --- | --- | --- |
| `f[i,s]` | cluster mean `CP` | Observed lineage DCF for clone `i`, sample `s` |
| `B[i,j]` | `B` | Clonal parent-child edge |
| `x[i,s]` | `x` | Exclusive clone prevalence |
| `y[i,s]` | `y` | Inferred lineage DCF, including descendants |
| `omega[i,site]` | `omega_site_vars` | Observed clone-site presence |
| `alpha[i,site]` | `alpha_vars` | Recursive admissibility |
| `mu[i,site]` | `mu_vars` | Clone-site pair selected in the refined tree |
| `d[i,j,site]` | `descent_edge_vars[child,parent,site]` | Descent edge at a site; implementation indices are child first |
| `M[i,source,destination]` | `migration_arc_node` | Migration edge for one clone in the refinement ILP |
| `G[source,destination]` | `location_arc` | Directed site pair used by at least one migration |
| `r` | `num_migrations` | `sum(mu[i,site], i > 0) - n` |
| `E` | `TREE_SCORE_NO_REWARD` | Weighted squared DCF fit error |

The default `--migration-presence-source node --node-trigger-mode proxy_all_obs` uses observed exclusive prevalence `f[i,s] - sum(B[i,j] * f[j,s])`, as in Supplementary Eq. 11. At least one supporting sample is needed per site; threshold indicators use the inherited numerical tolerance of `1e-10`.

Only **recursive admissibility** is implemented: `alpha[i,site]` holds precisely when the pair is observed or a child is admissible there, and `omega <= mu <= alpha`. The old direct-trigger admissibility mode and its standalone inference paths have been removed. The healthy root has no anatomical site or migration edges. `--max-root-out-degree 1` enforces one founding clone; `0` permits multiple root children.

Other presence signals remain available for sensitivity analyses. `--migration-presence-source lineage` uses inferred `y`; `observed` uses observed lineage DCF. Node trigger `any_x`/`avg_x` uses fitted exclusive `x` (site aggregation remains an OR across samples); `proxy` subtracts inferred child lineage DCFs. These choices change local presence evidence while retaining recursive admissibility. Use the defaults for the paper formulation.

## Main parameters

Run `python citupm.py --help` for all accepted options.

| Option | Default | Purpose |
| --- | --- | --- |
| `-f`, `--frequenciesFile` | required | Mutation/cluster DCF table |
| `-o`, `--pathOutputFilePrefix` | required | Output prefix; create its parent directory first |
| `--objective` | `cluster` | Cluster-mean fit or per-mutation fit (`mutation`) |
| `--clusterWeight` | `size` | `equal`, `size`, or `sqrt` cluster weights |
| `--sites-file` | omitted | Enable joint migration inference with sample/site mapping |
| `--presence-threshold` | `0.05` | Base presence threshold; paper protocol uses `0.025` |
| `--migration-presence-threshold` | base threshold | Override migration evidence threshold |
| `--max-root-out-degree` | `1` | Root degree cap; `0` is unlimited |
| `--metastatic-sites` | omitted | Sites excluded as root-child origins |
| `--minimize-migrations-at-score E_hat` | omitted | Minimize migrations under a DCF error bound |
| `--migration-score-epsilon` | `0.05` | Square-root error tolerance in Eq. 9 |
| `--migration-score-source` | `without-reward` | Bound residual error; `with-reward` is an alternative |
| `--migration-score-abs-tolerance` | omitted | Use `E <= reference + tolerance` instead of Eq. 9 |
| `--no-migration-score-tiebreak` | off | Skip stage 3; returned fit need not be best at that migration count |
| `--num-migrations`, `--max-migrations` | omitted | Require exactly this migration count |
| `--migration-weight` | `0` | Optional `gamma * r` objective penalty; must be zero in score-bounded mode |
| `--max-clone-site-pairs` | omitted | Cap `sum(mu)`; this is a pair count, not a migration count |
| `--numSolutions K` | omitted | Up to `K` best distinct clonal trees |
| `--enum-max-trees` | `20` | Near-optimal enumeration limit; `0` or `-1` means unlimited |
| `--enum-fit-abs`, `--enum-fit-rel` | `0`, `0` | Fit-window tolerances for near-optimal enumeration |
| `--max-refined-solutions` | `1` | Minimum-comigration refinements per clonal tree; `0` means all |
| `--fit-gap`, `--fit-gap-abs` | `0`, `0` | Gurobi relative/absolute fit optimality gaps |
| `--fit-time-limit` | omitted | Fit budget in seconds; shared by migration minimization and stage 3 |
| `--refinement-time-limit` | omitted | Limit in seconds per refinement solve |
| `--gurobi-max-time` | omitted | Default time limit per model |
| `--gurobi-threads` | solver default | Solver threads |
| `--gurobi-presolve` | `1` | Conservative presolve; `-1`, `0`, `2` select automatic/off/aggressive |

`--max-refined-solution` remains a singular alias. Additional options retain ancestry guidance (`-A`, `--autoAncestryTolerance`, `--disable-auto-ancestry-caps`), cluster mapping (`--fixClusterMapping`, on by default; `--no-fixClusterMapping`), filtering, clipping, optional rewards, fixed inputs, warm starts, and model diagnostics (`--write-mps`, `--write-mps-only`, `--debug-iis`).

`--min-mutations-per-cluster k` removes clusters with **at most** `k` mutations; default `0` retains nonempty clusters. `--filter-unbiased-var-ge v` removes clusters whose sum of within-sample unbiased variances is at least `v`. `--legacy-fit` retains the earlier fit objective solely for reproducing historical results; it does not select another admissibility model and is not used in the paper workflow.

An ancestry matrix is an unlabelled `(n+1) x (n+1)` numeric table in node order, including root `0`: `1` requires ancestry, `0` forbids it, and `-1` leaves it unspecified. Apply the same ancestry restrictions in reference and joint runs.

## Outputs

All files start with the `-o` prefix. Each `tree_N.txt` begins with alternating key/value lines, followed by blank-line-separated tables.

| File suffix | Contents |
| --- | --- |
| `.tree_N.txt` | Fit scores, counts/status, clonal edges, cluster/mutation assignments, `x`/`y`, mutation residuals, selected clone-site pairs and refined edges |
| `.tree_N_refined_K.txt` | Additional refinement `K` for topology `N` |
| `.tree_N_best_migration.txt` | Convenience copy selected by migration count, comigration count, then fit |
| `.treeScores.tsv` | One row per clonal topology; `numMigrations`, `numComigrations`, fit, refinement status and best-solution flags |
| `.cluster_variance.tsv` | Cluster sizes, variance statistics and filtering decisions |
| `.refined_nodes_solN.tsv` | `node, site, presence_signal, observed, mu, d, is_origin` |
| `.refined_edges_solN.tsv` | `from_node, from_site, to_node, to_site, type`; `type` is `descent` or `migration` |
| `.refined_solutions_solN.tsv` | Selected `mu` and origin indicators for each emitted refinement |
| `.refined_solution_summary_solN.tsv` | Per-refinement clone-site pair count, migration count, comigration count, site arcs and DAG flag |
| `.location_graphs_solN.tsv` | Directed site arcs for each emitted refinement |
| `.migration_score_optimization.tsv` | Score-bound optimization audit: reference, tolerance, exact cap, achieved score, selected count and proof status |
| `.core.gurobi.log`, `.refined_enum_solN.gurobi.log` | Solver diagnostics |

Sidecars for additional refinements and the convenience copy have matching `_refined_K` and `_best_migration` suffixes. Clone-site tables exclude the healthy root. Root descent edges have an empty source-site field; all such edges refer to the same healthy root. `presence_signal` is the maximum local evidence across samples at the site, rather than a fitted lineage DCF. With the default trigger it is observed exclusive prevalence and can be negative. `mu` records the selected presence; `d` and `is_origin` identify the clone's origin site.

`NUM_MIGRATIONS` is the number of migration edges, `sum(mu) - n`. `NUM_COMIGRATIONS` is the number of distinct directed source/destination site pairs in the optimal acyclic refinement, corresponding to the paper's comigration count in the absence of reseeding. For the bundled example these counts are **2** and **1**. If no optimal acyclic refinement is obtained, the program reports the integrated assignment with a non-`ok` `MIGRATION_STATUS` and `NUM_COMIGRATIONS = NA`; inspect the status and DAG flag before interpreting it.

`REFINED_ENUMERATION_STATUS = capped` with `REFINED_ENUMERATION_EXHAUSTIVE = NO` is expected when the requested refinement limit is reached, including the default limit of one. It does not mean that the returned minimum-comigration solution is suboptimal. `complete` means all minimum-comigration refinements were exhausted.

Scores are written at full precision so they can be reused as reference values. The raw Python/Gurobi solution-object dump has been removed; normal solver logs and a short completion message remain.

The four redundant `MIGRATION_SCORE_*`/`MIGRATION_BOUNDED_SCORE` values have been removed from tree headers. Their meanings were: selected score type, reference error `E_hat`, allowed upper error `(sqrt(E_hat)+epsilon)^2`, and achieved bounded error. Their audit information remains in `.migration_score_optimization.tsv`, outside the primary tree report. For example, reference `0.0004` with epsilon `0.05` gives an upper error of `0.0049`; achieved error `0.0004172275` satisfies that bound. Assignment indices/counts/deduplication bookkeeping and the duplicate `SUM_M_MINUS_I` metric are removed. Migration and comigration counts are reported directly. There are no personal-name prefixes in the new CLI, model variables, filenames or file sections; use the included plotters with these outputs.

## Plotting

Clonal tree PDFs, with sample/site legends and migration transitions:

```bash
python plot_with_gene_names.py results/patient/joint \
  --sample-sites /path/to/sites.tsv --with-site-transitions --with-node-info
```

Optionally add `--gene-info /path/to/annotations.tsv` for gene symbols, genomic locations and cancer roles (see `--help`), or `--site-colors /path/to/colors.tsv` with columns `site color`. `--output` sets a PDF path in single-file mode; `--output-dir` sets a destination for batch plots. A single `.refined_edges_solN.tsv` can also be passed to draw its migration graph; `--refined-nodes` and `--refined-tree` override inferred companion paths.

Expanded clone-site migration graphs, site-level comigration graphs, and conjoined summaries:

```bash
python plot_comigration_graphs.py results/patient/joint --sites-file /path/to/sites.tsv
```

This writes `.migration_graph.dot/pdf`, `.comigration_graph.dot/pdf/tsv`, and `.comigration_conjoined_solN.dot/pdf/tsv`. Graph edges are labelled by the number of migrating clones. Dotted node outlines mark selected pairs with `observed=0, mu=1`, i.e. unobserved seeding clones. `--migration-style detailed` includes clone/site text. `--no-conjoined` skips the conjoined plot; it sums counts over alternative refinements and must not be interpreted as one inferred history. By default, both plotters skip `_best_migration` copies to avoid duplicate figures; add `--include-best-migration-copy` to include them.

## Validation

```bash
python -m unittest discover -s tests -v
```

The integration checks use synthetic data to validate the fit/migration stages, recursive presence propagation, fixed migration counts, refinement enumeration, output schemas and plotting. A working Gurobi installation is required for solver checks; Graphviz and plotting dependencies are required for PDF checks.
