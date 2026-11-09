# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.19.3
#   kernelspec:
#     display_name: Python (base)
#     language: python
#     name: base
# ---

import pandas as pd
import numpy as np

pareto_optima_file = "results/min_trees/tracerx/optima.tsv"

us = pd.read_csv(
    pareto_optima_file,
    names=["person", "migrations", "comigrations", "quad", "tree"],
    sep="\t",
)

them = pd.read_csv(
    "results/fixed_trees/tracerx/optima.tsv",
    names=["person", "migrations", "comigrations", "quad", "tree"],
    sep="\t",
)

both = us.merge(them, on="person", suffixes=("_us", "_them"))

mach2 = pd.read_csv(
    "results/summary/mach2_migrations.tsv",
    sep="\t",
    names=["person", "migrations_mach2", "comigrations_mach2", "is_dag"],
)

metient = pd.read_csv(
    "results/summary/metient_migrations.tsv",
    sep="\t",
    names=["person", "migrations_metient", "comigrations_metient"],
)

cmp = both.merge(mach2, on="person", how="left").merge(metient, on="person", how="left").set_index("person", drop=False)


def migrations_cmp(a, b):
    if a < b:
        return "migrations_better"
    elif a > b:
        return "migrations_worse"
    else:
        return "migrations_same"


def mach2_comigrations_cmp(a, b):
    if a < b:
        return "0mach2_worse_cm"
    elif a > b:
        return "2mach2_better_cm"
    else:
        return "1mach2_same_cm"


def tree_score_cmp(a, b):
    if a < b:
        return "score_better"
    elif a > b:
        return "score_worse"
    else:
        return "score_same"


def tree_score_relcmp(a, b):
    if a < 0.95 * b:
        return "score_better"
    elif a > 1.05 * b:
        return "score_worse"
    else:
        return "score_same"


def mach2_cmp(a, b):
    if a < b:
        return "0mach2_worse"
    elif a > b:
        return "2mach2_better"
    else:
        return "1mach2_same"


def metient_cmp(a, b):
    if a < b:
        return "0metient_worse"
    elif a > b:
        return "2metient_better"
    else:
        return "1metient_same"


cmp["migrations_cmp"] = cmp.apply(
    lambda r: migrations_cmp(r["migrations_us"], r["migrations_them"]), axis=1
)

cmp["fixed_mig_cmp"] = cmp.apply(
    lambda r: migrations_cmp(r["migrations_them"], r["migrations_mach2"]),
    axis=1,
)

cmp["fixed_metient_cmp"] = cmp.apply(
    lambda r: migrations_cmp(r["migrations_them"], r["migrations_metient"]),
    axis=1,
)

cmp["mach2_cmp"] = cmp.apply(
    lambda r: mach2_cmp(r["migrations_us"], r["migrations_mach2"]), axis=1
)

cmp["metient_cmp"] = cmp.apply(
    lambda r: metient_cmp(r["migrations_us"], r["migrations_metient"]), axis=1
)

cmp["mach2_metient_cmp"] = cmp.apply(
    lambda r: metient_cmp(r["migrations_mach2"], r["migrations_metient"]), axis=1
)

cmp["mach2_comigrations_cmp"] = cmp.apply(
    lambda r: mach2_comigrations_cmp(
        r["comigrations_us"], r["comigrations_mach2"]
    ),
    axis=1,
)

cmp["tree_score_cmp"] = cmp.apply(
    lambda r: tree_score_cmp(r["quad_us"], r["quad_them"]), axis=1
)

cmp["tree_score_relcmp"] = cmp.apply(
    lambda r: tree_score_relcmp(r["quad_us"], r["quad_them"]), axis=1
)

# ## Absolute score comparisons

print(pd.crosstab(cmp["migrations_cmp"], cmp["tree_score_cmp"]).to_markdown())

print(pd.crosstab(cmp["mach2_cmp"], cmp["tree_score_cmp"]).to_markdown())

print(pd.crosstab(cmp["metient_cmp"], cmp["tree_score_cmp"]).to_markdown())

# ## Relative Score comparisions

print(
    pd.crosstab(cmp["migrations_cmp"], cmp["tree_score_relcmp"]).to_markdown()
)

print(pd.crosstab(cmp["mach2_cmp"], cmp["tree_score_relcmp"]).to_markdown())

print(pd.crosstab(cmp["metient_cmp"], cmp["tree_score_relcmp"]).to_markdown())

# ## There are no cases on the **fixed** CONIPHER tree that MACH2 or Metient gets fewer migrations than CITUP

cmp.fixed_mig_cmp.value_counts()

cmp.fixed_metient_cmp.value_counts()

# ## Cases where we have worse tree score

print(
    cmp.query('tree_score_cmp == "score_worse"')[
        [
            "person",
            "migrations_us",
            "quad_us",
            "migrations_them",
            "quad_them",
            "migrations_mach2",
        ]
    ].to_markdown(index=False)
)

# ## There are no cases in the joint pipeline where MACH2 or Metient gets fewer migrations

cmp.mach2_cmp.value_counts()

cmp.metient_cmp.value_counts()

# ## Metient is usually worse than MACH2, but sometimes gets fewer migrations

cmp.mach2_metient_cmp.value_counts()

print(
    cmp.query('mach2_cmp == "0mach2_worse" & tree_score_cmp == "score_same"')[
        ["person", "migrations_us", "migrations_mach2"]
    ].to_markdown(index=False)
)

# ## Mach2 often does not produce DAGS

mach2.is_dag.value_counts()

# ## Sometimes we don't produce DAGs

print(cmp.query('comigrations_us.isna()')[['person', 'migrations_us', 'quad_us', 'migrations_mach2', 'is_dag']].to_markdown(index=False))

# ## Our commigrations may be worse

# +

print(
    cmp.query(
        'mach2_comigrations_cmp == "2mach2_better_cm" & mach2_cmp == "1mach2_same" &tree_score_cmp == "score_same"'
    )[
        [
            "person",
            "migrations_us",
            "comigrations_us",
            "comigrations_mach2",
            "is_dag"
        ]
    ].to_markdown(index=False)
)
# -

# ## On the CONIPHER tree, how often do we produce DAGs

conipher_optima_dag = np.where(~cmp.comigrations_them.isna(), 'CITUP_DAG', 'CITUP_notDAG')

mach2_optima_dag = np.where(cmp.is_dag, 'MACH2_DAG', 'MACH2_notDAG')

print(pd.crosstab(conipher_optima_dag, mach2_optima_dag).to_markdown())

print(
    cmp.query(
        'comigrations_mach2 < comigrations_them & migrations_mach2 == migrations_them'
    )[
        [
            "person",
            "migrations_us",
            "comigrations_us",
            "comigrations_mach2",
            "is_dag"
        ]
    ].to_markdown(index=False)
)

print(
    cmp.query(
        'comigrations_mach2 < comigrations_them & migrations_mach2 == migrations_them'
    )[
        [
            "person",
            "migrations_them",
            "comigrations_them",
            "migrations_mach2",
            "comigrations_mach2",
            "is_dag"
        ]
    ].to_markdown(index=False)
)

print(
    cmp.query(
        'comigrations_them.isna()'
    )[
        [
            "person",
            "migrations_them",
            "migrations_mach2",
            "comigrations_mach2",
            "is_dag"
        ]
    ].to_markdown(index=False)
)

print(cmp.query('quad_us < quad_them & migrations_us < migrations_them')[['person', 'quad_us', 'quad_them', 'migrations_us', 'migrations_them', 'migrations_mach2', 'migrations_metient', 'is_dag']].to_markdown())

print(cmp.query('quad_us < quad_them')[['person', 'quad_us', 'quad_them', 'migrations_us', 'migrations_them', 'migrations_mach2', 'migrations_metient', 'is_dag']].to_markdown())

print(cmp.query('migrations_us != migrations_them & quad_us < quad_them')[['person', 'quad_us', 'quad_them', 'migrations_us', 'migrations_them', 'migrations_mach2', 'migrations_metient', 'is_dag']].to_markdown())

ccmp = cmp.copy()

ccmp['diff'] = ccmp.quad_them / (ccmp.quad_us + 1e-6)

ccmp.sort_values('diff', ascending=False)[['person', 'diff', 'quad_them', 'quad_us', 'migrations_us', 'migrations_them']].iloc[0:15,:]


