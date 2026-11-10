# CITUPm

CITUPm jointly infers a clonal tree, clone prevalences, and a migration history from mutation-cluster descendant cell fractions (DCFs) measured across tumor samples.  The nodes of a clonal tree, except the healthy root, represent tumor clones, and each of them corresponds to a distinct cluster of mutations.  A migration history is a tree that includes migrations between tumor sites, where there may be one or more samples from a site.

This repository contains the current inference entry point, `citupm.py`, and the two main plotting scripts. It also contains a pipeline for processing multiple cases in a cohort using CITUPm; see [CITUPm pipeline](CITUPm_pipeline/README.md).

## Setup

To download CITUPm and install the required Python dependencies in a Python virtual environment, execute the following code.

```bash
git clone https://github.com/algo-cancer/CITUPm.git
cd CITUPm
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

The virtual environment `CITUPm/.venv` needs to be active when you use CITUPm.  Once activated within a shell, it will stay active until you exit.  Thereafter, you should reactivate it using the following code.

```bash
cd CITUPm
source .venv/bin/activate
```

Install the Graphviz application separately so `dot` is on `PATH`; the Python `graphviz` package alone does not supply it. For example, use `sudo apt-get install graphviz` on Debian/Ubuntu, `brew install graphviz` on macOS, or your cluster's Graphviz module. See [Graphviz installation](https://graphviz.org/download/).

CITUPm was implemented using Gurobi optimization software, and a Gurobi license is required. Configure your institution's license or set `GRB_LICENSE_FILE` to your license file. See [Gurobi Python installation and licensing](https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python).

Once installation is complete, the following code should execute without error and optimize an example dataset.

```bash
python -c 'import gurobipy, pandas, graphviz; print(gurobipy.gurobi.version())'
dot -V
python citupm.py --help
bash examples/run_paper_workflow.sh
```
The example writes solver outputs and tree plots to `results/example/`.  If you wish, pass a different output directory as the script's first argument.

## Inputs

### Mutation frequencies (`-f`, required)

A whitespace-delimited table with a header. The first column contains unique mutation IDs, and the last column contains cluster IDs. Each intervening column contains the DCF values for the mutations in one sample.  Include only tumor-associated somatic mutations.

```text
mutID   primary_1   primary_2   metastasis_1   clusterID
m1      0.80        0.75        0.70           C1
m2      0.82        0.77        0.72           C1
m3      0.40        0.35        0.30           C2
m4      0.42        0.37        0.32           C2
```

DCFs are proportions, not percentages. Use purity-adjusted or unadjusted values consistently across the whole matrix. CITUPm does not infer purity or perform copy-number correction. `--clip-observed-at 1.0` changes values greater than 1.0 to 1.0 before fitting.

For each sample, the mean DCF of the mutations in a cluster supplies that cluster’s observed DCF.  With `--clusterWeight equal`, every cluster-sample pair has equal weight in the computation of DCF fit error, as defined in the supplied CITUPm RECOMB 2027 manuscript. With the script default `size`, a cluster's weight is its mutation count.

### Sample locations (`--sites-file`)

To activate migration history computation, supply the `--sites-file <file>` option with a sites file as described in this section. The file should be a two-column sample-to-site mapping, with an optional `sample site` header. Every sample named in the column headers of the frequency-table file must appear in the sites file. Multiple samples may map to one site; a cluster is considered present at a site if it is present in any sample from that site.  Presence in a sample is defined in the supplied CITUPm RECOMB 2027 manuscript. Sample and site labels must not contain whitespace.

```text
sample          site
primary_1       primary
primary_2       primary
metastasis_1    metastasis
```

`--metastatic-sites metastasis` optionally prevents that site from being the origin of the tumor (or tumors if more than one site of origin is permitted, which is only the case if `--max-root-out-degree` is not 1). Replace `metastasis` with one or more site names, separated by commas and/or spaces. Names must occur in the mapping, and at least one site must remain eligible as an origin.

### Fixed trees and warm starts

You may supply a fixed tree to infer an optimal DCF fit, or supply a tree as a starting point for migration minimization.  Trees are specified by two files.  The first file, a 'B' file, lists parent-child edges using integer node IDs.   Node `0` is the healthy root.  The second file, a 'delta' file, maps each non-root node from 'B' to a cluster name.

- `--fixed-B-file`: complete `parent child` edge list using node IDs; optional header `parent child`.
- `--fixed-delta-file`: map of each non-root node to a distinct cluster.  The first column contains cluster names, and the second contains node IDs; optional header `cluster node`.
- `--warm-start-first-phase-b` and `--warm-start-first-phase-delta`: the same formats, but used as a warm start rather than a fixed tree.  The warm start is only used when the `--minimize-migrations-at-score` option is passed.

## Paper workflow

The terminology and equations here follow the supplied RECOMB 2027 manuscript. Its default procedure has three optimization stages:

1. Minimize DCF fit error without migration constraints. Call the minimum DCF fit error `E_hat`.
2. Minimize migration count `r` subject to `E_r <= (sqrt(E_hat) + epsilon)^2`, where `E_r` is the minimal DCF fit error for `r` migrations (or infinity if there is no feasible tree with `r` migrations).  Here, `epsilon` is the tolerance specified by `--migration-score-epsilon`.  Let `r_star` be the smallest migration count satisfying this bound.
3. Fix the selected migration count `r_star` and minimize DCF fit error again. Stages 2 and 3 run in a single invocation by default.

After fitting, CITUPm seeks migration histories represented by trees that are consistent with the clonal tree and connect tumor sites.  We call these trees refined solutions.  We seek refined solutions that don't involve reseeding.  That is, no site may seed itself through a sequence of migrations.  Among those that do not involve reseeding, we seek to minimize the number of comigrations. For migration trees without reseeding the comigration count is the number of distinct directed source/destination site pairs in the refined tree. `--max-refined-solutions` controls how many such refined trees may be emitted. A value of `0` requests enumeration of all qualifying minimum-comigration refinements without reseeding for each clonal tree.

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
Stages 2 and 3 share the `--fit-time-limit` budget.

The [CITUPm pipeline](CITUPm_pipeline/README.md) implements the same three-stage procedure but for multiple cases in a cohort.



### Other inference modes

Optimize DCF fit error for an **exact** migration count, for example to explore the fit/migration tradeoff for a dataset:

```bash
python citupm.py -f /path/to/frequencies.tsv -o results/fixed_r/patient \
  --sites-file /path/to/sites.tsv --num-migrations 5 \
  --objective cluster --clusterWeight equal --clip-observed-at 1.0 \
  --presence-threshold 0.025 --numSolutions 1
```

Create `results/fixed_r/` first. `--max-migrations` is accepted as the manuscript's spelling of this option; despite that spelling, the inherited formulation fixes the count by equality. `--num-migrations` makes that behavior explicit. It cannot be combined with `--minimize-migrations-at-score`.

`--numSolutions K` enumerates up to `K` distinct topologies ranked by DCF fit error, including ties. It may return fewer when there are no more feasible topologies or a resource limit intervenes.

## Main parameters

Run `python citupm.py --help` for all accepted options.

| Option                                 | Default        | Purpose                                                                                                                 |
| -------------------------------------- | -------------- | ----------------------------------------------------------------------------------------------------------------------- |
| `-f`, `--frequenciesFile`              | required       | Mutation/cluster DCF table                                                                                              |
| `-o`, `--pathOutputFilePrefix`         | required       | Output prefix; create its parent directory first                                                                        |
| `--objective`                          | `cluster`      | Cluster-mean fit or per-mutation fit (`mutation`)                                                                       |
| `--clusterWeight`                      | `size`         | `equal`, `size`, or `sqrt` cluster weights                                                                              |
| `--sites-file`                         | omitted        | Enable joint migration inference with sample/site mapping                                                               |
| `--presence-threshold`                 | `0.05`         | Base presence threshold as defined in equation (11) of the supplied RECOMB 2027 manuscript; the manuscript uses `0.025` |
| `--max-root-out-degree`                | `1`            | Root degree cap; `0` is unlimited                                                                                       |
| `--metastatic-sites`                   | omitted        | Sites excluded as root-child origins                                                                                    |
| `--minimize-migrations-at-score E_hat` | omitted        | Minimize migrations under a DCF error bound                                                                             |
| `--migration-score-epsilon`            | `0.05`         | Square-root error tolerance                                                                                             |
| `--num-migrations`, `--max-migrations` | omitted        | Require exactly this migration count                                                                                    |
| `--numSolutions K`                     | omitted        | Up to `K` best distinct clonal trees                                                                                    |
| `--max-refined-solutions`              | `1`            | Minimum-comigration refinements per clonal tree; `0` means all                                                          |
| `--fit-time-limit`                     | omitted        | Fit budget in seconds; shared by migration minimization and stage 3                                                     |
| `--clip-observed-at`                   | omitted        | Caps observed DCF values before fitting                                                                                 |
| `--fixed-B-file`                       | omitted        | Supplies a fixed tree’s parent–child edges                                                                              |
| `--fixed-delta-file`                   | omitted        | Supplies the cluster-to-node mapping for a fixed tree                                                                   |
| `--warm-start-first-phase-b`           | omitted        | Supplies parent–child edges for a warm start                                                                            |
| `--warm-start-first-phase-delta`       | omitted        | Supplies the cluster-to-node mapping for a warm start                                                                   |
| `--gurobi-threads`                     | Solver default | Sets the number of solver threads                                                                                       |
| `--min-mutations-per-cluster`          | 0              | Removes clusters with at most the specified number of mutations                                                         |

`--min-mutations-per-cluster k` removes clusters with **at most** `k` mutations. The default value of `0` retains nonempty clusters.

## Outputs

### Filenames

All files start with the `-o` prefix.  Output files of interest to users are as follows.

| File suffix                                         | Contents                                                                                                                                                                      |
| --------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `.tree_N.txt`                                       | Fit scores, counts/status, clonal edges, cluster/mutation assignments, inferred node and lineage prevalences, mutation residuals, selected clone-site pairs and refined edges |
| `.tree_N_refined_K.txt`                             | Additional refinement `K` for topology `N`                                                                                                                                    |
| `.tree_N_best_migration.txt`                        | Convenience copy selected by migration count, comigration count, then fit                                                                                                     |
| `.treeScores.tsv`                                   | One row per clonal topology; `numMigrations`, `numComigrations`, fit, refinement status and best-solution flags                                                               |
| `.cluster_variance.tsv`                             | Cluster sizes, variance statistics and filtering decisions                                                                                                                    |
| `.refined_edges_solN.tsv`                           | `from_node, from_site, to_node, to_site, type`; `type` is `descent` or `migration`                                                                                            |
| `.location_graphs_solN.tsv`                         | Directed site arcs for each emitted refinement                                                                                                                                |
| `.core.gurobi.log`, `.refined_enum_solN.gurobi.log` | Solver diagnostics                                                                                                                                                            |

Filenames for additional refinements and the convenience copy have matching `_refined_K` and `_best_migration` suffixes. Clone-site tables exclude the healthy root. Root descent edges have an empty source-site field; all such edges refer to the same healthy root.

### Tree files (.tree_N.txt)

Each `tree_N.txt` begins with alternating key/value lines, followed by blank-line-separated tables.

`NUM_MIGRATIONS` is the number of migration edges. `NUM_COMIGRATIONS` is the number of distinct directed source/destination site pairs in an optimal refinement without reseeding, corresponding to the paper's comigration count. For the bundled example these counts are **2** and **1**. If no optimal refinement without reseeding is found, the program reports a non-`ok` `MIGRATION_STATUS` and `NUM_COMIGRATIONS = NA`. It emits an arbitrary refined tree that provides a migration history consistent with the clonal tree and has a minimal number of migrations for that clonal tree, provided the clonal tree is feasible.

`REFINED_ENUMERATION_STATUS = capped` with `REFINED_ENUMERATION_EXHAUSTIVE = NO` is expected when the requested refinement limit is reached, including the default limit of one. It does not mean that the returned minimum-comigration solution is suboptimal. `complete` means all minimum-comigration refinements were exhausted.

If optimality has not been proved and a feasible tree has been found, CITUPm will return the best feasible tree found so far.   To check optimality, check `FIT_OPTIMAL` (stage 1 is optimal), `MIGRATION_MINIMUM_PROVEN` (stage 2 is optimal), and `MIGRATION_SCORE_TIEBREAK_PROVEN` (stage 3 is optimal) for the value `YES`.

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

This writes `.migration_graph.dot/pdf`, `.comigration_graph.dot/pdf/tsv`, and `.comigration_conjoined_solN.dot/pdf/tsv`. Graph edges are labelled by the number of migrating clones. Dotted node outlines mark clones whose presence at a site has been inferred but not observed. The conjoined plot sums counts over alternative refinements and must not be interpreted as one inferred history.

## Validation

```bash
python -m unittest discover -s tests -v
```

The integration checks use synthetic data to validate the fit/migration stages, fixed migration counts, refinement enumeration, output schemas and plotting. A working Gurobi installation is required for solver checks; Graphviz and plotting dependencies are required for PDF checks.
