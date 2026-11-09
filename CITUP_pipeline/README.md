
## Tree optimization and migration analysis

This project uses Snakemake to run CITUPm tree optimization and migration
analyses for each person in a selected cohort. The workflow supports fitting
supplied fixed trees, optimizing trees using the no-migrations configuration,
minimizing migrations using a stage 1 reference score, and scanning specified
migration limits.

### Workflow targets

The aggregate targets run the corresponding analysis across the cohort.

| Target                      | Purpose                                                                           |
| --------------------------- | --------------------------------------------------------------------------------- |
| `stage0_fixed_trees_all`    | Fit mutation data using supplied fixed-tree inputs.                               |
| `stage1_optimal_trees_all`  | Optimize trees using the cohort’s no-migrations configuration.                    |
| `stage23_minimal_trees_all` | Minimize migrations using a reference score and warm starts derived from stage 1. |
| `migration_scan_all`        | Run separate optimizations for the migration limits listed in the samplesheet.    |

**The default target is** `**stage23_minimal_trees_all**`**.** It depends on
stage 1 results, so Snakemake schedules the required stage 1 jobs automatically.
Stage 0 and migration scans are independent branches and must be requested
separately.

### Requirements and input files

The workflow requires Snakemake, Python, pandas, and CITUPm optimizer
with its Gurobi dependency.

Run commands from the project root, with the following files available:

| Path                                     |Required for                               |
| ---                                      | ---                                       |
| `config.yaml`                            | Workflow configuration                    |
| `samplesheets/{cohort}/samplesheet.tsv`  | Cohort membership and per-person settings |
| `data/{cohort}/{person}/mutations.tsv`   | All processing rules                      |
| `data/{cohort}/{person}/fixed_B.tsv`.    | Stage 0                                   |
| `data/{cohort}/{person}/fixed_delta.tsv` | Stage 0                                   |
| `data/{cohort}/{person}/sites.tsv`       | Stages 2–3 and migration scans            |

The Snakefile does not generate these source data files or define their column
formats. Prepare them in the formats required by the configured optimizer.

### Samplesheet

The samplesheet is tab-separated and uses the following columns:

| Column       | Description|
| ---          | --- |
| `person`     | Person identifier matching the corresponding data directory. Use one row per person.|
| `tests`      | Comma-separated integer migration limits, such as `0,1,2,3`. Leave blank to omit that person from the aggregate migration scan.|
| `metastases` | Optional value passed to `--metastatic-sites`. Leave blank to omit this option; otherwise use the optimizer’s expected syntax.|

Include the column headers even when optional values are blank. The `tests`
column is parsed when the workflow loads, including when running targets other
than migration scans. The samplesheet reader treats `#` as a comment marker.

Each nonblank `tests` entry creates one scan job per listed limit. Limits are
paired with their corresponding person rather than combined across the cohort.

### Configuration

Set `cohort` in `config.yaml` to select the samplesheet and data directories.

Optimizer settings are organized under three sections, each keyed by cohort name:

|Configuration section. | Used by        |
|---                    | ---            |
|`citupm_no_migrations` | Stages 0 and 1 |
|`citupm_min_migrations`| Stages 2–3     |
|`citupm_scan`          | Migration scans|

Each applicable cohort entry must define:

- `citup_py`: path to the optimizer’s Python script.
    
- `args`: additional command-line arguments.
    
- `num_solutions`: value passed to `--numSolutions`.
    
- `dag_solutions`: value passed to `--max-refined-solutions` in stage 0, stages 2–3, and migration scans.
    

Stage 1 reads `dag_solutions` from the configuration but does not pass it to the optimizer.

### Running the workflow

Preview the default workflow and its shell commands:

```
snakemake --cores 8 --dry-run --printshellcmds
```

Run the default stages 1–3 analysis:

```
snakemake --cores 8 stage23_minimal_trees_all
```

Run another analysis by selecting its target:

```
# Fit supplied fixed trees.
snakemake --cores 8 stage0_fixed_trees_all

# Run stage 1 only.
snakemake --cores 8 stage1_optimal_trees_all

# Run the migration limits specified in the samplesheet.
snakemake --cores 8 migration_scan_all
```

### Outputs and logs

| Analysis       | Output directory|
| ---            | ---                                                          |
| Stage 0        | `results/fixed_trees/{cohort}/{person}/trees/`               |
| Stage 1        | `results/optimal_trees/{cohort}/{person}/trees/`             |
| Stages 2–3     | `results/min_trees/{cohort}/{person}/trees/`                 |
| Migration scan | `results/migration_scan/{cohort}/mig_{migrations}/{person}/` |

The aggregate scan target pads migrations to at least three digits, for example `mig_002`.

For stages 0, 1, and 2–3, logs are stored in `logs/citup.log` alongside each
person’s `trees` directory. Migration-scan logs are stored at:

```
results/migration_scan/{cohort}/logs/mig_{migrations}/{person}/citup.log
```

Each processing rule records the SLURM job ID, or `NA` when unavailable, and
redirects optimizer output and errors to its log.

### Computing resources

Every processing rule declares 8 threads and `mem_mb=8192`. The thread
allocation is passed to the optimizer through `--gurobi-threads`.

| Analysis       | Requested runtime   |
| ---            | ---                 |
| Stage 0        | 4 hours             |
| Stage 1        | 24 hours            |
| Stages 2–3     | 24 hours 30 minutes |
| Migration scan | 24 hours            |

These are resource declarations; their enforcement depends on the execution
backend. For cluster execution, configure an appropriate Snakemake executor or
profile. Recording a SLURM job ID does not itself submit jobs to SLURM.
