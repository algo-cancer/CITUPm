
# CITUPm pipeline

This project uses Snakemake to run CITUPm tree optimization and migration
analyses for each person in a selected cohort. The workflow supports fitting
supplied fixed trees, optimizing trees using the no-migrations configuration,
and minimizing migrations subject to a DCF fit error tolerance around the
solution from the no-migrations configuration.  

As a pipeline, it is designed for bulk processing rather than processing an
individual sample.  It takes as input a samplesheet specifying persons in a
cohort and a structured set of input files and directories.   It produces a
structured set of output files and directories within the directory `results/`.

As a Snakemake pipeline, it can be run in many differing compute environments,
from a single machine to nodes on a cluster.  Typically, Snakemake is invoked
in a Linux or Unix-like environment like MacOS. It is most practical to run this
pipeline on a cluster. We ran on the heterogenous NIH Biowulf cluster with a
Slurm scheduler using their Snakemake instructions.  We do not have any specific
instructions for how to run on your cluster; Snakemake documents the general
process but we expect that your cluster will have specific instructions.

Installation of CITUPm should include all specific software needed to run this
pipeline.   In general, the workflow requires Snakemake, Python, pandas, and
CITUPm optimizer with its Gurobi dependency.


### Workflow targets

Snakemake uses aggregate targets to specify what analyses should be performed.
The aggregate targets run the corresponding analysis across the cohort, which is
by default 60 simulations adapted from https://github.com/raphael-group/machina.

| Target                      | Purpose                                                                           |
| --------------------------- | --------------------------------------------------------------------------------- |
| `stage0_fixed_trees_all`    | Fit mutation data using supplied fixed-tree inputs.                               |
| `stage1_optimal_trees_all`  | Optimize trees using the cohort’s no-migrations configuration.                    |
| `stage23_minimal_trees_all` | Minimize migrations using a reference score and warm starts derived from stage 1. |

**The default target is** `**stage23_minimal_trees_all**`**.** Which is the complete 
CITUPm three stage pipeline (stage 1 is run as a dependancy.)

Run commands from the project root, with the following files available:

| Path                                     |Required for                               |
| ---                                      | ---                                       |
| `config.yaml`                            | Workflow configuration                    |
| `samplesheets/{cohort}/samplesheet.tsv`  | Cohort membership and per-person settings |
| `data/{cohort}/{person}/mutations.tsv`   | All processing rules                      |
| `data/{cohort}/{person}/fixed_B.tsv`.    | Stage 0                                   |
| `data/{cohort}/{person}/fixed_delta.tsv` | Stage 0                                   |
| `data/{cohort}/{person}/sites.tsv`       | Stages 2–3           |

The Snakefile does not generate these source data files or define their column
formats. Prepare them in the formats required by CITUPm.

### Samplesheet

The samplesheet is tab-separated and uses the following columns:

| Column       | Description|
| ---          | --- |
| `person`     | Person identifier matching the corresponding data directory. Use one row per person.|
| `metastases` | Optional value passed to `--metastatic-sites`. Leave blank to omit this option; otherwise use the optimizer’s expected syntax.|

### Configuration

Set `cohort` in `config.yaml` to select the cohort to run the pipeline on.  This
configuration value then determines the samplesheet, the paths to the the input
files as described above, and the paths to the output files as described below. 

Optimizer settings are organized under two sections, each keyed by cohort name:

|Configuration section. | Used by        |
|---                    | ---            |
|`citupm_no_migrations` | Stages 0 and 1 |
|`citupm_min_migrations`| Stages 2–3     |

Each applicable cohort entry must define:

- `citup_py`: path to the optimizer’s Python script.
    
- `args`: additional command-line arguments.
    
- `num_solutions`: the number of clonal trees T for CITUPm to compute.
    
- `dag_solutions`: the number of refined trees to be produced for each clonal
  tree T.  Set to zero to enumerage all refined trees that do not involve
  reseeding.
    

### Running the workflow

As discussed above, running the workflow is most practical on a cluster.  For
instance, there are 60 simulations in the default cohort and these would take a
considerable amount of time to run sequentially.  Nevertheless, we cannot know
the exact commmand that you would need to use on your cluster, so we give the
correct command lines to run sequentially.

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
profile.
