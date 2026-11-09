#!/usr/bin/env python3
"""CITUPm clonal-tree and migration inference. See README.md for the paper workflow."""
from gurobipy import *
import argparse, math, os, sys, time

epsilon = 1e-10

def roundInt(x): return int(x + 0.5)
def isFloat(v):
    try: float(v); return True
    except: return False
def floatToStr(x, d=4): return ("{:." + str(d) + "f}").format(float(x)) if isFloat(x) else x

def safe_val(var):
    """Return variable value even if a later stage failed (uses X or Xn)."""
    for attr in ('X', 'Xn'):
        try:
            return var.getAttr(attr)
        except Exception:
            continue
    return 0.0

def roundVar(var):
    return roundInt(safe_val(var))

def is_dag_directed(edge_list):
    """Return True iff directed graph given by edge_list has no cycles."""
    if not edge_list:
        return True
    nodes = set()
    adj = {}
    for u, v in edge_list:
        nodes.add(u); nodes.add(v)
        if u not in adj:
            adj[u] = []
        adj[u].append(v)
    color = {n: 0 for n in nodes}  # 0=unvisited, 1=visiting, 2=done

    def dfs(n):
        color[n] = 1
        for nb in adj.get(n, []):
            c = color.get(nb, 0)
            if c == 1:
                return False
            if c == 0 and not dfs(nb):
                return False
        color[n] = 2
        return True

    for n in nodes:
        if color[n] == 0 and not dfs(n):
            return False
    return True

def migration_graph_dag_flag(sol):
    """Return 'YES'/'NO'/'NA' depending on whether the migration graph is a DAG."""
    edges = sol.get("refined_edges")
    if edges is None and sol.get("num_migrations") is None:
        return "NA"
    if not edges:
        return "YES"
    mig_edges = []
    for (fn, fs, tn, ts, etype) in edges:
        if etype == "migration":
            mig_edges.append((fs, ts))
    return "YES" if is_dag_directed(mig_edges) else "NO"


# ---------------- CLI ----------------
parser = argparse.ArgumentParser(
    description='CITUPm: joint clonal tree and migration history inference',
    add_help=True
)

# Required
parser.add_argument('-f','--frequenciesFile', required=True)
parser.add_argument('-o','--pathOutputFilePrefix', required=True)

# Optional ancestry
parser.add_argument('-A','--ancestryMatrixFile', type=str, default=None)
parser.add_argument('--autoAncestryTolerance', type=float, default=None)

# Objective / data options
parser.add_argument('--objective', choices=['mutation','cluster'], default='cluster')
parser.add_argument('--clusterWeight', choices=['size','equal','sqrt'], default='size')
parser.add_argument('--presence-threshold', type=float, default=0.05)
parser.add_argument('--co-presence-reward', type=float, default=0.0)
parser.add_argument('--duo-present-reward', type=float, default=0.0)
parser.add_argument('--duo-absent-reward', type=float, default=0.0)
parser.add_argument('--duo-presence-source', choices=['observed', 'inferred', 'node', 'stacking'], default='observed',
                    help='Presence signal used for duo rewards: observed (omega), inferred (y), node (x), or stacking (CITUPm site-level mu calls).')
parser.add_argument('--duo-presence-threshold', type=float, default=None,
                    help='Threshold for duo presence when source is inferred/node (defaults to --presence-threshold).')
parser.add_argument('--min-mutations-per-cluster', type=int, default=0)
parser.add_argument('--filter-unbiased-var-ge', type=float, default=None,
                    help='Filter OUT clusters whose perSampleVarSum_unbiased >= threshold (across-sample sum of unbiased within-sample variances).')
mapping_group = parser.add_mutually_exclusive_group()
mapping_group.add_argument('--fixClusterMapping', action='store_true', default=True)
mapping_group.add_argument('--no-fixClusterMapping', dest='fixClusterMapping', action='store_false',
                    help='Optimize the cluster-node mapping instead of fixing the default mapping. '
                         'An explicit --fixed-delta-file still fixes the mapping.')

# Gurobi (global + per-stage)
parser.add_argument('--gurobi-threads', type=int, default=None)
parser.add_argument('--gurobi-max-time', type=float, default=None)
parser.add_argument('--gurobi-presolve', type=int, choices=[-1, 0, 1, 2], default=1,
                    help='Gurobi presolve level: -1=automatic, 0=off, 1=conservative, '
                         '2=aggressive (default: 1).')
parser.add_argument('--write-mps', type=str, default=None)
parser.add_argument('--write-mps-only', action='store_true', default=False,
                    help='Write the stage-1 model requested by --write-mps and exit before optimization. '
                         'Intended for reproducible solver diagnostics.')

# Stage-1 solver controls (fit)
parser.add_argument('--fit-gap', type=float, default=0.0)
parser.add_argument('--fit-gap-abs', type=float, default=0.0)
parser.add_argument('--fit-time-limit', type=float, default=None)

# Optional migration minimization with a user-supplied tree-score reference.
parser.add_argument('--minimize-migrations-at-score', type=float, default=None,
                    help='Minimize migrations subject to a tree-score bound around this reference. '
                         'By default the score excludes rewards and the bound follows paper Eq. (9).')
parser.add_argument('--migration-score-source', choices=['with-reward', 'without-reward'],
                    default='without-reward',
                    help='Score bounded in migration minimization: residual only '
                         '(without-reward, default), or residual minus rewards (with-reward).')
score_tiebreak_group = parser.add_mutually_exclusive_group()
score_tiebreak_group.add_argument('--migration-score-tiebreak', dest='migration_score_tiebreak',
                    action='store_true', default=True,
                    help='After minimizing migrations, fix that count and minimize the selected '
                         'tree score (default; applies to --minimize-migrations-at-score).')
score_tiebreak_group.add_argument('--no-migration-score-tiebreak', dest='migration_score_tiebreak',
                    action='store_false',
                    help='Skip the internal tree-score minimization after minimizing migrations. '
                         'Keep the migration incumbent and continue refinement within the score bound.')
score_tolerance_group = parser.add_mutually_exclusive_group()
score_tolerance_group.add_argument('--migration-score-epsilon', type=float, default=0.05,
                    help='Eq. (9) tolerance: score <= (sqrt(reference) + epsilon)^2 '
                         '(default: 0.05; requires a nonnegative reference).')
score_tolerance_group.add_argument('--migration-score-abs-tolerance', type=float, default=None,
                    help='Use an additive score bound, score <= reference + tolerance, instead '
                         'of Eq. (9). Supports negative reward-adjusted reference scores.')

# Enumeration over near-optimal trees
parser.add_argument('--enum-fit-rel', type=float, default=0.0,
                    help='Enumerate trees with fit ≤ best*(1+rel) (or ≤ best+abs).')
parser.add_argument('--enum-fit-abs', type=float, default=0.0,
                    help='Enumerate trees with fit ≤ best+abs (combined with rel).')
parser.add_argument('--enum-max-trees', type=int, default=20,
                    help='Max number of trees to collect (including the best). Set 0 or -1 for unlimited (subject to --enum-time-limit).')
parser.add_argument('--enum-time-limit', type=float, default=None,
                    help='Overall time cap (seconds) for tree enumeration.')
# Number of distinct clonal trees to enumerate (top-K by fit).
# When set, enumerate successive best distinct topologies; ties are allowed.
parser.add_argument('--numSolutions', dest='num_solutions', type=int, default=None,
                    help='Enumerate up to K distinct clonal trees ranked by DCF fit error. Overrides enum-max-trees and near-optimal fit bounds; score-bounded mode retains its cap and migration count.')
parser.add_argument('--sites-file', type=str, default=None,
                    help='Sample↦site mapping to enable migration constraints in the core ILP.')
parser.add_argument('--metastatic-sites', nargs='*', default=None,
                    help='Optional metastatic site names, separated by commas and/or spaces. '
                         'These sites cannot be selected as root-child descent origins '
                         '(paper healthy-root model: the root itself has no site).')
parser.add_argument('--migration-weight', type=float, default=0.0,
                    help='Weight gamma on migration count r = sum(mu[v,site]) - n in the core objective.')
parser.add_argument('--max-clone-site-pairs', type=int, default=None,
                    help='Optional cap on selected clone-site pairs: sum(mu[v,site]) <= M.')
parser.add_argument('--num-migrations', '--max-migrations', dest='num_migrations', type=int, default=None,
                    help='Exact migration count (also accepted as manuscript --max-migrations): sum_{non-root v,site} mu[v,site] - n == value.')
parser.add_argument('--migration-presence-threshold', type=float, default=None,
                    help='Threshold on chosen CITUPm presence signal (defaults to --presence-threshold).')
parser.add_argument('--migration-presence-source', choices=['observed', 'node', 'lineage'], default='node',
                    help='Presence evidence: observed lineage DCF, node signal, or inferred lineage DCF (y). The default node/proxy_all_obs combination implements paper Eq. (11).')
parser.add_argument('--node-trigger-mode', choices=['any_x', 'avg_x', 'proxy', 'proxy_all_obs'], default='proxy_all_obs',
                    help='How node-source builds per-sample presence before site OR: any_x/avg_x use x, proxy uses CP*delta-child_mass via By, proxy_all_obs uses CP*delta-child_mass via observed child signal.')
parser.add_argument('--refinement-time-limit', type=float, default=None,
                    help='Optional time limit (seconds) for each CITUPm migration ILP solve.')
parser.add_argument('--max-refined-solutions', '--max-refined-solution',
                    dest='max_refined_solutions', type=int, default=1, metavar='N',
                    help='Maximum number of minimum-comigration refined trees to emit per clonal tree. '
                         'Default: 1. Set 0 for all minimum-comigration refinements.')
parser.add_argument('--max-root-out-degree', '--max-root-out-deg', dest='max_root_out_degree',
                    type=int, default=1,
                    help='Maximum allowed out-degree of healthy root node 0 (default: 1). '
                         'Set 0 for unlimited.')
parser.add_argument('--fixed-B-file', type=str, default=None,
                    help='Optional fixed tree topology as edge list file with columns: parent child (node indices).')
parser.add_argument('--fixed-delta-file', type=str, default=None,
                    help='Optional fixed cluster-node assignment file with columns: cluster node (non-root node index).')
parser.add_argument('--warm-start-first-phase-b', '--warm-start-fisrst-phase-b',
                    dest='warm_start_first_phase_b', type=str, default=None,
                    help='Optional B MIP start for the first --minimize-migrations-at-score solve. '
                         'Uses the --fixed-B-file format (parent child), but adds no constraints.')
parser.add_argument('--warm-start-first-phase-delta', '--warm-start-fisrst-phase-delta',
                    dest='warm_start_first_phase_delta', type=str, default=None,
                    help='Optional delta MIP start for the first --minimize-migrations-at-score solve. '
                         'Uses the --fixed-delta-file format (cluster node), but adds no constraints. '
                         'Use --no-fixClusterMapping if the mapping should be optimized.')

parser.add_argument('--disable-auto-ancestry-caps', action='store_true', default=False,
                    help='Ignore CP-based allowed-edge mask (A_allowed_all) when building B/q/migration variables.')

# Debug
parser.add_argument('--debug-iis', action='store_true', default=False)

# Optional preprocessing of observed frequencies
parser.add_argument('--clip-observed-at', type=float, default=None,
                    help='Clip observed mutation/cluster frequencies at this value before fitting '
                         '(applies to both objectives).')

# Reproduce legacy (buggy) fit behavior across both objectives
parser.add_argument('--legacy-fit', action='store_true', default=False,
                    help='Reproduce legacy mutation fit that used z on all nodes with M=1. '
                         'Applies this residual for BOTH objectives so results match the legacy output.')

# Mutation objective linearization control
parser.add_argument('--mut-abs-bigM', type=float, default=None,
                    help='Override absolute big-M for mutation objective linearization. '
                         'Default: auto = max(1.0, max observed frequency).')

args = parser.parse_args()

score_bounded_migrations = args.minimize_migrations_at_score is not None
score_reference = args.minimize_migrations_at_score
score_upper = None
if not args.migration_score_tiebreak and not score_bounded_migrations:
    parser.error('--no-migration-score-tiebreak requires --minimize-migrations-at-score.')
if (args.warm_start_first_phase_b is not None or args.warm_start_first_phase_delta is not None) and not score_bounded_migrations:
    parser.error('First-phase warm starts require --minimize-migrations-at-score.')
if score_bounded_migrations:
    if not math.isfinite(score_reference):
        parser.error('--minimize-migrations-at-score must be finite.')
    if args.sites_file is None:
        parser.error('--minimize-migrations-at-score requires --sites-file.')
    if args.num_migrations is not None:
        parser.error('--minimize-migrations-at-score cannot be combined with '
                     '--num-migrations, which fixes the migration count.')
    if not math.isfinite(args.migration_weight) or args.migration_weight != 0.0:
        parser.error('--minimize-migrations-at-score requires --migration-weight 0; '
                     'migration count is minimized separately from the tree score.')
    if args.migration_score_source == 'without-reward' and score_reference < 0:
        parser.error('A without-reward reference score cannot be negative.')
    if args.migration_score_abs_tolerance is not None:
        tolerance = args.migration_score_abs_tolerance
        if not math.isfinite(tolerance) or tolerance < 0:
            parser.error('--migration-score-abs-tolerance must be finite and nonnegative.')
        score_upper = score_reference + tolerance
    else:
        tolerance = args.migration_score_epsilon
        if not math.isfinite(tolerance) or tolerance < 0:
            parser.error('--migration-score-epsilon must be finite and nonnegative.')
        if score_reference < 0:
            parser.error('Eq. (9) is undefined for a negative reference score. '
                         'Use --migration-score-abs-tolerance for an additive bound, '
                         'or --migration-score-source without-reward with a nonnegative reference.')
        norm_upper = math.sqrt(score_reference) + tolerance
        score_upper = norm_upper * norm_upper
    if not math.isfinite(score_upper):
        parser.error('The migration score bound must be finite.')

def _log_args(parsed):
    """Display all CLI arguments so runs are traceable."""
    print("=== CITUPm arguments ===")
    for key in sorted(vars(parsed)):
        print(f"{key} = {getattr(parsed, key)}")
    print("=== end arguments ===")

_log_args(args)

def parse_optional_site_list(raw_values, option_name):
    if raw_values is None:
        return []
    parsed = []
    seen = set()
    for raw_value in raw_values:
        for site in str(raw_value).replace(",", " ").split():
            site = site.strip()
            if site and site not in seen:
                parsed.append(site)
                seen.add(site)
    if not parsed:
        raise ValueError(f"{option_name} was provided but no sites were listed.")
    return parsed

presence_threshold = max(0.0, min(1.0, float(args.presence_threshold)))
migration_presence_threshold = presence_threshold if args.migration_presence_threshold is None else max(0.0, float(args.migration_presence_threshold))
migration_presence_source = args.migration_presence_source
node_trigger_mode = args.node_trigger_mode
migration_weight = float(args.migration_weight)
max_clone_site_pairs_arg = None if args.max_clone_site_pairs is None else max(0, int(args.max_clone_site_pairs))
num_migrations_arg = None if args.num_migrations is None else max(0, int(args.num_migrations))
max_refined_solutions_arg = max(0, int(args.max_refined_solutions))
max_root_out_degree_arg = int(args.max_root_out_degree)
if max_root_out_degree_arg < 0:
    raise ValueError("--max-root-out-degree must be nonnegative (0 for unlimited).")
if max_root_out_degree_arg == 0:
    max_root_out_degree_arg = None
metastatic_sites_arg = parse_optional_site_list(args.metastatic_sites, "--metastatic-sites")
metastatic_sites = []
metastatic_site_set = set()
co_presence_reward = float(args.co_presence_reward)
duo_present_reward = float(args.duo_present_reward)
duo_absent_reward = float(args.duo_absent_reward)
duo_presence_source = args.duo_presence_source
duo_presence_threshold = presence_threshold if args.duo_presence_threshold is None else max(0.0, min(1.0, float(args.duo_presence_threshold)))
min_mut_per_cluster = max(0, int(args.min_mutations_per_cluster))

if duo_presence_source == 'stacking':
    if abs(duo_absent_reward) > epsilon:
        raise ValueError("--duo-presence-source stacking only supports --duo-present-reward; Section 6 does not define a duo-absent term.")

# ---------------- IO ----------------
f, c = {}, {}
mutIDs, clusterIDs = [], []
clusterID_counts = {}
with open(args.frequenciesFile, "r") as fh:
    hdr = fh.readline().strip().split()
    sampleIDs = [s.strip() for s in hdr[1:-1]]
    numSamples = len(sampleIDs)
    for line in fh:
        cols = line.strip().split()
        assert len(cols) - 2 == numSamples, "Bad row (sample count mismatch): " + line
        mid = cols[0]
        assert mid not in mutIDs, f"Duplicate mutation id: {mid}"
        mutIDs.append(mid)
        f[mid] = {sampleIDs[i]: float(cols[i+1]) for i in range(numSamples)}
        cid = cols[-1]
        if cid not in clusterIDs: clusterIDs.append(cid)
        c[mid] = cid
# Build cluster->mut list and snapshots
cluster_to_muts = {cid: [] for cid in clusterIDs}
for mid in mutIDs:
    cluster_to_muts[c[mid]].append(mid)

# Optional: clip observed frequencies before computing CP/weights, so both objectives see the same targets
if args.clip_observed_at is not None:
    cap = float(args.clip_observed_at)
    for mid in mutIDs:
        for s in sampleIDs:
            if f[mid][s] > cap:
                f[mid][s] = cap

# Big-M for mutation objective linearization: must dominate |f - y| when delta==0.
# y in [0,1]. If user provided --mut-abs-bigM, use it; otherwise auto = max(1.0, max observed f).
try:
    _f_max = max((f[mid][s] for mid in mutIDs for s in sampleIDs), default=1.0)
except Exception:
    _f_max = 1.0
OBS_DCF_UB = max(1.0, float(_f_max))
if args.mut_abs_bigM is not None:
    BIGM_MUT_ABS = float(args.mut_abs_bigM)
    if BIGM_MUT_ABS < _f_max - 1e-9:
        print(f"[warn] --mut-abs-bigM={BIGM_MUT_ABS} < max observed frequency {_f_max:.4f}; mutation objective may not match cluster objective.")
else:
    BIGM_MUT_ABS = OBS_DCF_UB

orig_map = dict(c)
orig_mut = list(mutIDs)
orig_clu = list(clusterIDs)

# Apply min-mutations-per-cluster (strict >) first
min_k = max(0, int(args.min_mutations_per_cluster))
kept_min = [cid for cid in orig_clu if len(cluster_to_muts[cid]) > min_k]

# Compute unbiased per-sample variance sum for each original cluster
def _unbiased_sum(cid):
    muts = cluster_to_muts[cid]
    k = len(muts)
    if k <= 1:
        return 0.0
    total_pop = 0.0
    for sname in sampleIDs:
        vals = [f[m][sname] for m in muts]
        mu = (sum(vals)/float(k)) if k else 0.0
        var_pop = (sum((v-mu)**2 for v in vals)/float(k)) if k else 0.0
        total_pop += var_pop
    return total_pop * (float(k)/float(k-1))

unb = {cid: _unbiased_sum(cid) for cid in orig_clu}

# If requested, filter by unbiased variance >= threshold
thr = args.filter_unbiased_var_ge
if thr is not None:
    thr = float(thr)
    retained_clusters = [cid for cid in kept_min if unb[cid] < thr]
else:
    retained_clusters = kept_min

if not retained_clusters:
    raise ValueError('All clusters filtered out by --min-mutations-per-cluster and/or --filter-unbiased-var-ge; consider lowering thresholds.')

# Emit unbiased variance report and concise stdout
report_path = args.pathOutputFilePrefix + '.cluster_variance.tsv'
with open(report_path, 'w') as fv:
    fv.write('clusterID	k	perSampleVarSum_unbiased	threshold	kept\n')
    thr_str = '' if args.filter_unbiased_var_ge is None else str(thr)
    for cid in orig_clu:
        k = len(cluster_to_muts[cid])
        kept = 1 if cid in retained_clusters else 0
        fv.write(str(cid)+'\t'+str(k)+'\t'+('%.6f' % unb[cid])+'\t'+thr_str+'\t'+str(kept)+'\n')
        print('Cluster '+str(cid)+': unbiasedVar=' + ('%.6f' % unb[cid]) + ' kept=' + str(kept))

# Restrict to retained clusters and rebuild counts/mappings
clusterIDs = retained_clusters
mutIDs = [mid for mid in orig_mut if orig_map[mid] in clusterIDs]
cluster_to_muts = {cid: [mid for mid in cluster_to_muts[cid] if mid in mutIDs] for cid in clusterIDs}
clusterID_counts = {cid: len(cluster_to_muts[cid]) for cid in clusterIDs}
c = {mid: orig_map[mid] for mid in mutIDs}

numMutations = len(mutIDs)
numClusters  = len(clusterIDs)
treeSize = numClusters + 1

# cluster prevalences (means) after filtering
CP = {cid: {s: (0.0 if len(cluster_to_muts[cid])==0 else
               sum(f[m][s] for m in cluster_to_muts[cid]) / float(len(cluster_to_muts[cid])))
           for s in sampleIDs}
      for cid in clusterIDs}

omega_cluster = {cid: {s: 1 if CP[cid][s] >= presence_threshold else 0 for s in sampleIDs}
                 for cid in clusterIDs}

# Identify samples that retain at least one tumour cluster after filtering
sample_requires_tumor = {s: any(omega_cluster[cid][s] == 1 for cid in clusterIDs) for s in sampleIDs}

# CITUPm migration inputs (sample↦site)
# If stacking/metastatic-site constraints are requested without an explicit site map,
# assume one sample per site.
migration_enabled = (
    args.sites_file is not None or
    args.metastatic_sites is not None or
    (duo_presence_source == 'stacking' and abs(duo_present_reward) > epsilon)
)
sample_to_site = {}
site_to_samples = {}
site_list = []
if args.sites_file is not None:
    with open(args.sites_file, "r") as sf:
        for line in sf:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            toks = line.split()
            if len(toks) < 2:
                continue
            if toks[0].lower() == "sample" and toks[1].lower() == "site":
                continue
            sample, site = toks[0], toks[1]
            sample_to_site[sample] = site
            if site not in site_to_samples:
                site_to_samples[site] = []
                site_list.append(site)
            site_to_samples[site].append(sample)
    missing = [s for s in sampleIDs if s not in sample_to_site]
    if missing:
        raise ValueError("sites-file missing samples: " + ", ".join(missing))
    if not site_list:
        raise ValueError("sites-file parsed but no site entries were found.")
elif migration_enabled:
    for sample in sampleIDs:
        sample_to_site[sample] = sample
        site_to_samples[sample] = [sample]
        site_list.append(sample)

if metastatic_sites_arg:
    known_sites = set(site_list)
    unknown_metastatic_sites = [site for site in metastatic_sites_arg if site not in known_sites]
    if unknown_metastatic_sites:
        raise ValueError(
            "--metastatic-sites contains unknown sites: {}. Known sites are: {}".format(
                ", ".join(unknown_metastatic_sites),
                ", ".join(site_list)
            )
        )
    root_origin_sites = [site for site in site_list if site not in set(metastatic_sites_arg)]
    if not root_origin_sites:
        raise ValueError("--metastatic-sites marks every known site as metastatic; no root origin site remains.")
    metastatic_sites = list(metastatic_sites_arg)
    metastatic_site_set = set(metastatic_sites)
    print("metastatic_sites_normalized = " + ",".join(metastatic_sites))

ordered_clusters = list(clusterIDs)
clusterIndex = {cid: i+1 for i,cid in enumerate(ordered_clusters)}  # root=0
# ancestry input
A_input = [[-1]*treeSize for _ in range(treeSize)]
if args.ancestryMatrixFile:
    A_input = []
    with open(args.ancestryMatrixFile,"r") as amf:
        for line in amf:
            line=line.strip()
            if not line: continue
            A_input.append([int(float(x)) for x in line.split()])
    assert len(A_input)==treeSize and all(len(r)==treeSize for r in A_input), "Bad ancestry matrix size"

# ancestry feasibility from CP (optional)
A_allowed_all = [[1]*treeSize for _ in range(treeSize)]
if args.autoAncestryTolerance is not None:
    tol = float(args.autoAncestryTolerance)
    def node_to_cluster(v): return None if v==0 else ordered_clusters[v-1]
    for v in range(treeSize):
        for u in range(treeSize):
            if v==u: A_allowed_all[v][u] = 1; continue
            if u==0: A_allowed_all[v][u] = 0; continue
            if v==0: A_allowed_all[v][u] = 1; continue
            cv, cu = node_to_cluster(v), node_to_cluster(u)
            ok = all(CP[cv][s] + tol >= CP[cu][s] for s in sampleIDs)
            A_allowed_all[v][u] = 1 if ok else 0

# with open(args.pathOutputFilePrefix + ".A_allowed_all_samples.tsv","w") as fa:
#     for v in range(treeSize):
#         fa.write("\t".join(str(A_allowed_all[v][u]) for u in range(treeSize))+"\n")

# Precompute which directed edges are allowed (used to sparsify B/q and migration vars)
allowed_edge = [[False]*treeSize for _ in range(treeSize)]
for v in range(treeSize):
    for u in range(treeSize):
        if v == u:
            allowed_edge[v][u] = False
            continue
        forb = False
        # External ancestry forbids direct edge
        a = A_input[v][u]
        if a in (0, 1) and a == 0:
            forb = True
        # CP-based monotonicity mask forbids direct edge (unless disabled)
        if (not args.disable_auto_ancestry_caps) and (A_allowed_all[v][u] == 0):
            forb = True
        allowed_edge[v][u] = (not forb)

# Convenience list of allowed directed edges (v,u) with v != u
ALLOWED_PAIRS = [(v,u) for v in range(treeSize) for u in range(treeSize)
                 if v != u and allowed_edge[v][u]]

# Optional fixed-input mode: user can provide a fixed topology (B) and/or
# fixed cluster-to-node assignment (delta), and the solver will infer the rest.
def _iter_clean_lines(path):
    with open(path, "r") as fh:
        for lno, raw in enumerate(fh, 1):
            line = raw.strip()
            if (not line) or line.startswith("#"):
                continue
            yield lno, line

def _parse_two_cols(path, first_name, second_name):
    rows = []
    header_seen = False
    for lno, line in _iter_clean_lines(path):
        cols = line.replace(",", " ").split()
        if len(cols) < 2:
            raise ValueError(f"{path}:{lno}: expected 2 columns ({first_name} {second_name})")
        # Allow one header row, e.g. "parent child" or "cluster node".
        if (not header_seen) and (cols[0].lower() == first_name.lower()) and (cols[1].lower() == second_name.lower()):
            header_seen = True
            continue
        rows.append((lno, cols[0], cols[1]))
    return rows

def _load_fixed_delta(path):
    rows = _parse_two_cols(path, "cluster", "node")
    if not rows:
        raise ValueError(f"{path}: no data rows found")
    out = {}
    used_nodes = set()
    for lno, cid, node_txt in rows:
        if cid not in clusterIDs:
            raise ValueError(f"{path}:{lno}: unknown cluster '{cid}'")
        try:
            node = int(node_txt)
        except Exception:
            raise ValueError(f"{path}:{lno}: node must be an integer, got '{node_txt}'")
        if node <= 0 or node >= treeSize:
            raise ValueError(f"{path}:{lno}: node must be in [1, {treeSize-1}], got {node}")
        if cid in out:
            raise ValueError(f"{path}:{lno}: duplicate cluster '{cid}'")
        if node in used_nodes:
            raise ValueError(f"{path}:{lno}: node {node} assigned to multiple clusters")
        out[cid] = node
        used_nodes.add(node)

    missing = [cid for cid in clusterIDs if cid not in out]
    extra_nodes = sorted(set(range(1, treeSize)) - used_nodes)
    if missing:
        raise ValueError(f"{path}: missing clusters in fixed delta: {','.join(map(str, missing))}")
    if extra_nodes:
        raise ValueError(f"{path}: unassigned non-root nodes in fixed delta: {','.join(map(str, extra_nodes))}")
    return out

def _load_fixed_B(path, check_model_constraints=True):
    rows = _parse_two_cols(path, "parent", "child")
    if not rows:
        raise ValueError(f"{path}: no data rows found")
    edges = set()
    indeg = [0 for _ in range(treeSize)]
    children = {v: [] for v in range(treeSize)}
    for lno, p_txt, c_txt in rows:
        try:
            p = int(p_txt); ch = int(c_txt)
        except Exception:
            raise ValueError(f"{path}:{lno}: parent/child must be integers, got '{p_txt} {c_txt}'")
        if p < 0 or p >= treeSize or ch < 0 or ch >= treeSize:
            raise ValueError(f"{path}:{lno}: node ids must be in [0, {treeSize-1}]")
        if p == ch:
            raise ValueError(f"{path}:{lno}: self-edge ({p},{ch}) is invalid")
        if (p, ch) in edges:
            raise ValueError(f"{path}:{lno}: duplicate edge ({p},{ch})")
        if check_model_constraints and not allowed_edge[p][ch]:
            raise ValueError(
                f"{path}:{lno}: edge ({p},{ch}) is forbidden by ancestry/allowed-edge filters"
            )
        edges.add((p, ch))
        indeg[ch] += 1
        children.setdefault(p, []).append(ch)

    if len(edges) != treeSize - 1:
        raise ValueError(f"{path}: expected exactly {treeSize-1} edges, got {len(edges)}")
    if indeg[0] != 0:
        raise ValueError(f"{path}: root node 0 must have indegree 0")
    for v in range(1, treeSize):
        if indeg[v] != 1:
            raise ValueError(f"{path}: node {v} must have indegree 1, got {indeg[v]}")

    visited = set()
    stack = [0]
    while stack:
        cur = stack.pop()
        if cur in visited:
            continue
        visited.add(cur)
        for ch in children.get(cur, []):
            stack.append(ch)
    if len(visited) != treeSize:
        raise ValueError(f"{path}: edges do not form a rooted tree spanning all nodes")
    if check_model_constraints and max_root_out_degree_arg is not None:
        root_out_degree = sum(1 for p, _ in edges if p == 0)
        if root_out_degree > max_root_out_degree_arg:
            raise ValueError(
                f"{path}: root out-degree {root_out_degree} exceeds "
                f"--max-root-out-degree={max_root_out_degree_arg}"
            )
    return edges

fixed_B_edges = None
if args.fixed_B_file is not None:
    fixed_B_edges = _load_fixed_B(args.fixed_B_file)
    print(f"[info] Loaded fixed B with {len(fixed_B_edges)} edges from {args.fixed_B_file}")

fixed_delta_assign = None
if args.fixed_delta_file is not None:
    fixed_delta_assign = _load_fixed_delta(args.fixed_delta_file)
    print(f"[info] Loaded fixed delta with {len(fixed_delta_assign)} assignments from {args.fixed_delta_file}")

warm_start_B_edges = None
if args.warm_start_first_phase_b is not None:
    # Validate the tree format, but let Gurobi reject starts violating model constraints.
    warm_start_B_edges = _load_fixed_B(args.warm_start_first_phase_b, check_model_constraints=False)
    print(f'[info] Loaded first-phase warm-start B from {args.warm_start_first_phase_b} (not fixed).')

warm_start_delta_assign = None
if args.warm_start_first_phase_delta is not None:
    warm_start_delta_assign = _load_fixed_delta(args.warm_start_first_phase_delta)
    print(f'[info] Loaded first-phase warm-start delta from {args.warm_start_first_phase_delta} (not fixed by this file).')
    active_fixed_delta = fixed_delta_assign if fixed_delta_assign is not None else (clusterIndex if args.fixClusterMapping else None)
    if active_fixed_delta is not None and warm_start_delta_assign != active_fixed_delta:
        print('[warn] Warm-start delta differs from the fixed mapping and may be rejected. '
              'Use --no-fixClusterMapping and omit --fixed-delta-file to optimize delta.', file=sys.stderr)

def apply_first_phase_warm_start(core, b_edges=None, delta_assign=None):
    """Seed only B/delta; leave frequencies, migrations and other variables undefined."""
    seeded = 0
    if b_edges is not None:
        for edge, var in core['B'].items():
            var.Start = int(edge in b_edges)
            seeded += 1
    if delta_assign is not None:
        for (cid, node), var in core['delta'].items():
            var.Start = int(delta_assign[cid] == node)
            seeded += 1
    if seeded:
        print(f'[info] Applied first-phase warm start to {seeded} B/delta variables; '
              'Gurobi will complete or reject the start, without adding constraints.')

# ---------------- utils ----------------
def setup_env(model, log_suffix=""):
    if args.gurobi_threads is not None:
        try: model.setParam('Threads', max(1,int(args.gurobi_threads)))
        except: pass
    # nodefile/mem (optional via env)
    try:
        nfs = os.environ.get('GRB_NODEFILE_START')
        if nfs: model.setParam('NodefileStart', float(nfs))
    except: pass
    try:
        nfdir = os.environ.get('GRB_NODEFILE_DIR') or os.environ.get('SLURM_TMPDIR') or os.environ.get('TMPDIR')
        if nfdir: model.setParam('NodefileDir', nfdir)
    except: pass
    try:
        mem_mb = os.environ.get('GRB_MEMLIMIT_MB')
        if mem_mb: model.setParam('MemLimit', float(mem_mb))
    except: pass
    if args.gurobi_max_time is not None:
        try: model.setParam('TimeLimit', float(args.gurobi_max_time))
        except: pass
    try: model.setParam('LogFile', args.pathOutputFilePrefix + f"{log_suffix}.gurobi.log")
    except: pass
    # Conservative presolve is the default; aggressive presolve caused numerical
    # failures on the CRUK0311 exact-migration models.
    try: model.setParam('Presolve', int(args.gurobi_presolve))
    except: pass

# -------- CITUPm migration helpers --------
def build_tree_maps(edge_set):
    parent = {}
    children = {v: [] for v in range(treeSize)}
    for (p, ch) in edge_set:
        parent[ch] = p
        children.setdefault(p, []).append(ch)
    return parent, children

def _get_stage1_y(sol):
    y_vals = sol.get("stage1_y")
    if y_vals:
        return y_vals
    core_i = sol.get("core")
    if core_i is None:
        return {}
    return {(v, s): safe_val(core_i["y"][v, s]) for v in range(treeSize) for s in sampleIDs}

def _get_stage1_x(sol):
    x_vals = sol.get("stage1_x")
    if x_vals:
        return x_vals
    core_i = sol.get("core")
    if core_i is None:
        return {}
    return {(v, s): safe_val(core_i["x"][v, s]) for v in range(treeSize) for s in sampleIDs}

def _get_stage1_omega(sol):
    core_i = sol.get("core")
    if core_i is None:
        return {}
    return {(v, s): roundInt(safe_val(core_i["omega"][v, s])) for v in range(treeSize) for s in sampleIDs}

def _get_delta_assign(sol):
    delta_assign = sol.get("stage1_delta")
    if delta_assign:
        return delta_assign
    delta_assign = {}
    core_i = sol.get("core")
    if core_i:
        delta_vars = core_i["delta"]
        for cid in clusterIDs:
            for v in range(treeSize):
                if roundVar(delta_vars[cid, v]) == 1:
                    delta_assign[cid] = v
                    break
    return delta_assign


def _children_from_solution(sol, parent_map=None):
    children = {v: [] for v in range(treeSize)}
    if parent_map is not None:
        for child, parent in parent_map.items():
            children.setdefault(parent, []).append(child)
        return children

    edge_set = sol.get("b_edges", set())
    if not edge_set:
        core_i = sol.get("core")
        if core_i is not None:
            try:
                edge_set = capture_B(core_i)
            except NameError:
                edge_set = set()
    for parent, child in edge_set:
        children.setdefault(parent, []).append(child)
    return children

def _observed_lineage_by_node(sol):
    node_to_cluster = {}
    for cid, node in _get_delta_assign(sol).items():
        if node is not None:
            node_to_cluster[node] = cid

    observed = {}
    for v in range(treeSize):
        cid = node_to_cluster.get(v)
        for s in sampleIDs:
            observed[v, s] = CP.get(cid, {}).get(s, 0.0) if cid is not None else 0.0
    return observed

def compute_sample_signal_values(sol, parent_map=None):
    if migration_presence_source == 'observed':
        return _get_stage1_omega(sol)
    if migration_presence_source == 'lineage':
        return _get_stage1_y(sol)
    if migration_presence_source != 'node':
        return _get_stage1_x(sol)

    if node_trigger_mode in ('any_x', 'avg_x'):
        return _get_stage1_x(sol)

    observed_lineage = _observed_lineage_by_node(sol)
    children = _children_from_solution(sol, parent_map=parent_map)
    inferred_lineage = _get_stage1_y(sol)

    sig_vals = {}
    for v in range(treeSize):
        for s in sampleIDs:
            # The root has no observed cluster lineage. In the core CITUPm model,
            # root admissibility is handled separately and the root is not an
            # observed migration pair.
            if v == 0:
                sig_vals[v, s] = 0.0
                continue

            node_obs = observed_lineage.get((v, s), 0.0)
            if node_trigger_mode == 'proxy_all_obs':
                child_mass = sum(observed_lineage.get((child, s), 0.0)
                                 for child in children.get(v, []))
            else:
                child_mass = sum(inferred_lineage.get((child, s), 0.0)
                                 for child in children.get(v, []))
            sig_vals[v, s] = node_obs - child_mass
    return sig_vals

def compute_site_signal(sol, site_to_samples_map, parent_map=None):
    sig_vals = compute_sample_signal_values(sol, parent_map=parent_map)
    site_signal = {v: {} for v in range(treeSize)}
    for v in range(treeSize):
        for site, samples in site_to_samples_map.items():
            vals = [sig_vals.get((v, s), 0.0) for s in samples]
            site_signal[v][site] = max(vals) if vals else 0.0
    return site_signal

def allow_metastatic_root_child_origin_sites(allowed, parent_map, sites, root_node=0):
    """Paper model: metastatic restrictions affect d[0,j,site], not A or mu."""
    return allowed


def build_mu_assignment_bounds(site_signal, parent_map, sites, root_node=0, core=None):
    # In recursive mode, reuse the exact omega/alpha solution from the fixed-tree core
    # model. This avoids rebuilding threshold decisions from floating-point values.
    if core is not None:
        omega_vars = core.get("omega_site_vars")
        admissibility_vars = core.get("alpha_vars")
        if omega_vars is not None and admissibility_vars is not None:
            allowed = {
                (v, site)
                for v in range(1, treeSize)
                for site in sites
                if roundVar(admissibility_vars[v, site]) >= 1
            }
            observed = {
                (v, site)
                for v in range(1, treeSize)
                for site in sites
                if roundVar(omega_vars[v, site]) == 1
            }
            allow_metastatic_root_child_origin_sites(allowed, parent_map, sites, root_node=root_node)
            return allowed, observed

    observed = set()
    for v in range(1, treeSize):
        for site in sites:
            if not site_to_samples.get(site, []):
                continue
            if site_signal.get(v, {}).get(site, 0.0) >= migration_presence_threshold - 1e-12:
                observed.add((v, site))


    children = {v: [] for v in range(treeSize)}
    for child, parent in parent_map.items():
        children.setdefault(parent, []).append(child)

    omega = set(observed)

    allowed = set()
    visited = set()
    postorder = []
    stack = [(root_node, 0)]
    while stack:
        node, state = stack.pop()
        if state == 0:
            if node in visited:
                continue
            visited.add(node)
            stack.append((node, 1))
            for child in children.get(node, []):
                stack.append((child, 0))
        else:
            postorder.append(node)
    for v in range(treeSize):
        if v not in visited:
            postorder.append(v)

    for v in postorder:
        if v == root_node:
            continue
        for site in sites:
            if (v, site) in omega or any((child, site) in allowed for child in children.get(v, [])):
                allowed.add((v, site))
    allow_metastatic_root_child_origin_sites(allowed, parent_map, sites, root_node=root_node)
    return allowed, observed


def nonroot_pair_count_from_mu(mu_sol, sites, root_node=0):
    return sum(int(mu_sol.get((v, site), 0)) for v in range(treeSize) if v != root_node for site in sites)

def migration_count_from_mu(mu_sol, sites, root_node=0):
    total_pairs = nonroot_pair_count_from_mu(mu_sol, sites, root_node=root_node)
    return max(total_pairs - (treeSize - 1), 0)

def migration_count_expr_from_mu_vars(mu_vars, sites, root_node=0):
    return quicksum(mu_vars.get((v, site), 0) for v in range(treeSize) for site in sites if v != root_node) - (treeSize - 1)

def site_migration_edges_from_labeling(mu_sol, origin_site, sites, root_support=None, root_node=0):
    mig_edges = []
    for v in range(treeSize):
        if v == root_node:
            continue
        origin = origin_site.get(v)
        if origin is None:
            continue
        for site in sites:
            if site != origin and mu_sol.get((v, site), 0) == 1:
                mig_edges.append((origin, site))
    return mig_edges

def mu_bit_for_labeling(v, site, mu_sol, origin_site, root_support=None, root_node=0):
    if v == root_node:
        return 0
    return mu_sol.get((v, site), 0)

def build_refined_migration_outputs(parent_map, sites, mu_sol, origin_site, site_signal, observed, root_support=None):
    d_sol = {(v, site): (1 if v != 0 and origin_site.get(v) == site else 0) for v in range(treeSize) for site in sites}

    refined_edges = []
    for v in range(1, treeSize):
        origin = origin_site.get(v)
        if origin is None:
            continue
        for site in sites:
            if mu_bit_for_labeling(v, site, mu_sol, origin_site, root_support=root_support) == 1 and site != origin:
                refined_edges.append((v, origin, v, site, "migration"))
    for child, parent in parent_map.items():
        origin = origin_site.get(child)
        if origin is None:
            continue
        if parent == 0:
            refined_edges.append((parent, "", child, origin, "descent"))
        elif mu_sol.get((parent, origin), 0) == 1:
            refined_edges.append((parent, origin, child, origin, "descent"))

    refined_nodes = []
    for v in range(1, treeSize):
        for site in sites:
            refined_nodes.append(dict(
                node=v,
                site=site,
                presence_signal=site_signal.get(v, {}).get(site, 0.0),
                observed=1 if (v, site) in observed else 0,
                mu=mu_bit_for_labeling(v, site, mu_sol, origin_site, root_support=root_support),
                d=d_sol.get((v, site), 0),
                is_origin=1 if origin_site.get(v) == site else 0
            ))

    return d_sol, refined_edges, refined_nodes

def build_core_migration_fallback(sol, parent_map, sites, site_signal, observed,
                                  target_num_migrations, enum_result, root_node=0):
    core_i = sol["core"]
    mu_vars = core_i.get("mu_vars") or {}
    d_edge_vars = core_i.get("descent_edge_vars")
    mu_sol = {(v, site): roundVar(mu_vars[v, site])
             for v in range(1, treeSize) for site in sites}
    origin_site = extract_origin_labeling(
        root_node, parent_map, sites, mu_sol, d_edge_vars=d_edge_vars
    )
    d_sol, refined_edges, refined_nodes = build_refined_migration_outputs(
        parent_map, sites, mu_sol, origin_site, site_signal, observed,
        root_support=None
    )
    location_edges = sorted(set(site_migration_edges_from_labeling(
        mu_sol, origin_site, sites, root_support=None, root_node=root_node
    )))
    dag_flag = "YES" if is_dag_directed(location_edges) else "NO"
    return dict(
        migration_status="no_minimum_comigration_solution_core_assignment_reported",
        num_migrations=target_num_migrations,
        num_clone_site_pairs=nonroot_pair_count_from_mu(mu_sol, sites, root_node=root_node),
        num_comigrations=None,
        refined_solution_count=0,
        refined_enumeration_status=enum_result.get("status", "unknown"),
        refined_enumeration_exhaustive=enum_result.get("exhaustive", False),
        migration_location_graph_edges=location_edges,
        mu=mu_sol,
        origin_indicator=d_sol,
        origin_site=origin_site,
        migration_root_support=None,
        refined_edges=refined_edges,
        refined_nodes=refined_nodes,
        site_signal=site_signal,
        migration_tiebreak_status="refinement_failed_core_assignment_reported",
        migration_dag_before_postprocess=dag_flag,
        migration_dag_postprocess_status="not_attempted_refinement_failed",
        migration_labeling_source="core_integrated_assignment_fallback_non_refined",
        migration_relabel_change_count=0,
        migration_dag_enum_status="superseded_by_refinement_failed",
        migration_dag_enum_objective=None,
        migration_dag_solutions=[],
        refined_solutions=[]
    )


def extract_origin_labeling(root_node, parent_map, sites, mu_sol, d_root_vars=None, d_edge_vars=None):
    origin_site = {}
    origin_site[root_node] = None

    for child, parent in parent_map.items():
        chosen = None
        if d_edge_vars is not None:
            for site in sites:
                dvar = d_edge_vars.get((child, parent, site))
                if dvar is not None and roundVar(dvar) == 1:
                    chosen = site
                    break
        if chosen is None:
            for site in sites:
                if parent == root_node and mu_sol.get((child, site), 0) == 1:
                    chosen = site
                    break
                if parent != root_node and mu_sol.get((child, site), 0) == 1 and mu_sol.get((parent, site), 0) == 1:
                    chosen = site
                    break
        origin_site[child] = chosen

    return origin_site


def refined_tree_signature(mu_sol, origin_site, sites, root_node=0):
    """Signature of the Supplementary Section 2.3 decision variables mu and d/delta."""
    enumerated_nodes = [v for v in range(treeSize) if v != root_node]
    return (
        tuple(int(mu_sol.get((v, site), 0)) for v in enumerated_nodes for site in sites),
        tuple(origin_site.get(v) for v in enumerated_nodes)
    )

def enumerate_min_comigration_refined_trees(sol_idx, sol, parent_map, sites, site_signal,
                                             observed, target_num_migrations, max_solutions=0):
    """Enumerate minimum-comigration DAG refined trees for fixed B and migration count.

    This is the reduced ILP from Supplementary Section 2.3 of the manuscript. mu_vars encode non-root
    mu variables, origin_vars encode each non-root node's descent/origin site (including
    root-to-child choices d[0,j,site]), migration_arc_node encodes M_ist for
    non-root nodes, and location_arc encodes G_st.
    """
    root_node = sol["core"].get("root", 0)
    max_solutions = max(0, int(max_solutions))
    admissible, observed_pairs = build_mu_assignment_bounds(
        site_signal, parent_map, sites, root_node=root_node, core=sol["core"]
    )
    observed_pairs = set(observed_pairs)
    root_children = sorted(child for child, parent in parent_map.items()
                           if parent == root_node)

    model = Model(f"CITUPmRefinedEnum_{sol_idx}")
    setup_env(model, f".refined_enum_sol{sol_idx}")
    if args.refinement_time_limit is not None:
        try: model.setParam('TimeLimit', float(args.refinement_time_limit))
        except: pass

    mu_vars = {}
    origin_vars = {}
    for v in range(1, treeSize):
        for site in sites:
            allowed = (v, site) in admissible
            lb = 1 if (v, site) in observed_pairs else 0
            ub = 1 if allowed else 0
            if lb and not ub:
                return dict(status=f"infeasible_observed_pair_{v}_{site}", solutions=[],
                            objective=None, exhaustive=True)
            mu_vars[v, site] = model.addVar(vtype=GRB.BINARY, lb=lb, ub=ub,
                                            name=f"refined_mu[{v},{site}]")
            origin_vars[v, site] = model.addVar(vtype=GRB.BINARY, lb=0, ub=ub,
                                            name=f"refined_origin[{v},{site}]")
    model.update()

    for v in range(1, treeSize):
        model.addConstr(quicksum(origin_vars[v, site] for site in sites) == 1,
                        name=f"refined_origin_onehot[{v}]")
        for site in sites:
            model.addConstr(origin_vars[v, site] <= mu_vars[v, site],
                            name=f"refined_origin_subset[{v},{site}]")
    for child, parent in parent_map.items():
        for site in sites:
            if parent == root_node:
                if site in metastatic_site_set:
                    model.addConstr(origin_vars[child, site] == 0,
                                    name=f"refined_root_child_not_metastatic_origin[{child},{site}]")
            else:
                model.addConstr(
                    origin_vars[child, site] <= mu_vars[parent, site],
                    name=f"refined_descent_exists[{child},{site}]"
                )

    migration_count_expr = migration_count_expr_from_mu_vars(mu_vars, sites, root_node=root_node)
    model.addConstr(migration_count_expr == int(target_num_migrations),
                    name="refined_fixed_migration_count")
    if max_clone_site_pairs_arg is not None:
        model.addConstr(quicksum(mu_vars.values()) <= max_clone_site_pairs_arg,
                        name="refined_total_pairs_cap")

    # Stacking rewards can change when refinement chooses a different mu assignment.
    refined_stacking_expr = None
    fixed_reward = None
    if score_bounded_migrations and duo_presence_source == 'stacking' and abs(duo_present_reward) > epsilon:
        refined_stacking_expr = LinExpr()
        for child, parent in parent_map.items():
            if parent == root_node:
                continue
            for site in sites:
                overlap = model.addVar(vtype=GRB.BINARY,
                                       name=f"refined_stacking[{parent},{child},{site}]")
                model.addConstr(overlap <= mu_vars[parent, site])
                model.addConstr(overlap <= mu_vars[child, site])
                model.addConstr(overlap >= mu_vars[parent, site] + mu_vars[child, site] - 1)
                refined_stacking_expr += overlap
            refined_stacking_expr -= 1
        fixed_reward = sol['reward_stage1'] - duo_present_reward * sol['core']['duo_present_expr'].getValue()
        if args.migration_score_source == 'with-reward':
            refined_score_expr = sol['residual_stage1'] - fixed_reward - duo_present_reward * refined_stacking_expr
            model.addConstr(refined_score_expr <= score_upper, name='refined_migration_score_cap')

    migration_arc_node = {}
    location_arc = {}
    topo = {}
    if len(sites) >= 2:
        for site in sites:
            topo[site] = model.addVar(vtype=GRB.INTEGER, lb=0, ub=len(sites) - 1,
                                      name=f"refined_topo[{site}]")
        migration_nodes = [v for v in range(treeSize) if v != root_node]
        for v in migration_nodes:
            for source in sites:
                for destination in sites:
                    if source == destination:
                        continue
                    arc = model.addVar(vtype=GRB.BINARY,
                                       name=f"refined_Marc[{v},{source},{destination}]")
                    migration_arc_node[v, source, destination] = arc
                    model.addConstr(arc <= origin_vars[v, source],
                                    name=f"refined_Marc_le_D[{v},{source},{destination}]")
                    model.addConstr(arc <= mu_vars[v, destination],
                                    name=f"refined_Marc_le_M[{v},{source},{destination}]")
                    model.addConstr(arc >= origin_vars[v, source] + mu_vars[v, destination] - 1,
                                    name=f"refined_Marc_ge_and[{v},{source},{destination}]")
        for source in sites:
            for destination in sites:
                if source == destination:
                    continue
                gvar = model.addVar(vtype=GRB.BINARY,
                                    name=f"refined_G[{source},{destination}]")
                location_arc[source, destination] = gvar
                arc_terms = [migration_arc_node[v, source, destination]
                             for v in migration_nodes]
                for arc_index, arc in enumerate(arc_terms):
                    model.addConstr(gvar >= arc,
                                    name=f"refined_G_or_lb[{source},{destination},{arc_index}]")
                model.addConstr(gvar <= quicksum(arc_terms),
                                name=f"refined_G_or_ub[{source},{destination}]")
                model.addConstr(
                    topo[destination] >= topo[source] + 1 - len(sites) * (1 - gvar),
                    name=f"refined_location_dag[{source},{destination}]"
                )

    comigration_expr = quicksum(location_arc.values())
    model.setObjective(comigration_expr, GRB.MINIMIZE)
    model.Params.MIPGap = 0
    model.Params.MIPGapAbs = 0
    model.optimize()

    if model.Status != GRB.OPTIMAL or model.SolCount == 0:
        return dict(status=f"minimum_not_proven_status_{model.Status}", solutions=[],
                    objective=None, exhaustive=False, model=model)

    minimum_comigrations = roundInt(comigration_expr.getValue())
    model.addConstr(comigration_expr == minimum_comigrations,
                    name="refined_fix_minimum_comigrations")
    model.setObjective(0.0, GRB.MINIMIZE)
    model.update()

    refined_solutions = []
    seen = set()
    enum_status = "not_completed"
    exhaustive = False
    enumerated_mu_vars = [
        mu_vars[v, site]
        for v in range(1, treeSize)
        for site in sites
    ]
    enumerated_d_vars = [
        origin_vars[v, site]
        for v in range(1, treeSize)
        for site in sites
    ]
    all_decision_vars = enumerated_mu_vars + enumerated_d_vars
    while True:
        if max_solutions > 0 and len(refined_solutions) >= max_solutions:
            enum_status = "capped"
            exhaustive = False
            break
        model.optimize()
        solve_status = model.Status
        if solve_status == GRB.INFEASIBLE:
            enum_status = "complete"
            exhaustive = True
            break
        if solve_status not in (GRB.OPTIMAL, GRB.TIME_LIMIT):
            enum_status = f"stopped_status_{solve_status}"
            exhaustive = False
            break
        if model.SolCount == 0:
            enum_status = ("time_limit_no_incumbent" if solve_status == GRB.TIME_LIMIT
                           else f"no_incumbent_status_{solve_status}")
            exhaustive = False
            break

        mu_sol = {(v, site): roundVar(mu_vars[v, site])
                 for v in range(1, treeSize) for site in sites}
        origin_site = {}
        origin_site[root_node] = None
        for v in range(1, treeSize):
            origin_site[v] = next((site for site in sites if roundVar(origin_vars[v, site]) == 1), None)
        signature = refined_tree_signature(mu_sol, origin_site, sites, root_node=root_node)
        if signature in seen:
            enum_status = "duplicate_signature"
            exhaustive = False
            break
        seen.add(signature)

        d_sol, refined_edges, refined_nodes = build_refined_migration_outputs(
            parent_map, sites, mu_sol, origin_site, site_signal, observed_pairs,
            root_support=None
        )
        location_edges = sorted((source, destination)
                                for (source, destination), gvar in location_arc.items()
                                if roundVar(gvar) == 1)
        refined_solutions.append(dict(
            migration_status="ok",
            refined_solution_index=len(refined_solutions),
            num_migrations=migration_count_from_mu(mu_sol, sites, root_node=root_node),
            num_clone_site_pairs=nonroot_pair_count_from_mu(mu_sol, sites, root_node=root_node),
            num_comigrations=minimum_comigrations,
            migration_location_graph_edges=location_edges,
            mu=mu_sol,
            origin_indicator=d_sol,
            origin_site=origin_site,
            migration_root_support=None,
            refined_edges=refined_edges,
            refined_nodes=refined_nodes,
            site_signal=site_signal,
            migration_tiebreak_status="refinement_min_comigrations",
            migration_dag_before_postprocess="YES",
            migration_dag_postprocess_status="not_needed",
            migration_labeling_source="refinement_joint_enumeration",
            migration_relabel_change_count=0,
            migration_dag_enum_status="superseded_by_refinement",
            migration_dag_enum_objective=None,
            migration_dag_solutions=[]
        ))
        if refined_stacking_expr is not None:
            refined_reward = fixed_reward + duo_present_reward * roundInt(refined_stacking_expr.getValue())
            refined_solutions[-1]['reward_stage1'] = refined_reward
            refined_solutions[-1]['fit'] = sol['residual_stage1'] - refined_reward

        selected = [var for var in all_decision_vars if roundVar(var) == 1]
        unselected = [var for var in all_decision_vars if roundVar(var) == 0]
        model.addConstr(quicksum(1 - var for var in selected) +
                        quicksum(var for var in unselected) >= 1,
                        name=f"refined_nogood[{len(refined_solutions)}]")
        model.update()
        # Adding the exclusion cut can reset the model's optimization status.
        if solve_status == GRB.TIME_LIMIT:
            enum_status = "time_limit_partial"
            exhaustive = False
            break

    for refined in refined_solutions:
        refined["refined_solution_count"] = len(refined_solutions)
        refined["refined_enumeration_status"] = enum_status
        refined["refined_enumeration_exhaustive"] = exhaustive
    return dict(status=enum_status, solutions=refined_solutions,
                objective=minimum_comigrations, exhaustive=exhaustive, model=model)

# --------------- build joint core model ---------------
def build_core():
    root=0
    m=Model('CITUPm_core'); setup_env(m,".core")
    # delta
    delta={}
    for cid in clusterIDs:
        for v in range(treeSize):
            delta[cid,v]=m.addVar(vtype=GRB.BINARY, name=f"delta[{cid},{v}]")
    # each cluster ↦ exactly one node
    for cid in clusterIDs:
        m.addConstr(quicksum(delta[cid,v] for v in range(treeSize))==1)
    # node occupancy: no cluster at root; exactly one per non-root
    for v in range(treeSize):
        if v==root:
            for cid in clusterIDs: m.addConstr(delta[cid, v]==0)
        else:
            m.addConstr(quicksum(delta[cid,v] for cid in clusterIDs)==1)
    # optionally fix mapping
    if fixed_delta_assign is not None:
        for cid in clusterIDs:
            vfix = fixed_delta_assign[cid]
            m.addConstr(delta[cid, vfix]==1)
            for v in range(treeSize):
                if v!=vfix: m.addConstr(delta[cid,v]==0)
    elif args.fixClusterMapping:
        for cid in ordered_clusters:
            vfix = clusterIndex[cid]
            m.addConstr(delta[cid, vfix]==1)
            for v in range(treeSize):
                if v!=vfix: m.addConstr(delta[cid,v]==0)

    # omega[v,s]
    omega={}
    for v in range(treeSize):
        for s in sampleIDs:
            omega[v,s]=m.addVar(vtype=GRB.BINARY, name=f"omega[{v},{s}]")
    # root absent (presence indicator)
    for s in sampleIDs: m.addConstr(omega[root,s]==0)
    # non-root omega = presence of its assigned cluster in that sample
    for v in range(1,treeSize):
        for s in sampleIDs:
            contrib=[delta[cid,v] for cid in clusterIDs if omega_cluster[cid][s]==1]
            if contrib: m.addConstr(omega[v,s]==quicksum(contrib))
            else:       m.addConstr(omega[v,s]==0)
    # no sample ONLY root present: require tumour cluster only when at least one survives filtering
    for s in sampleIDs:
        if sample_requires_tumor.get(s, False):
            m.addConstr(quicksum(omega[v,s] for v in range(1,treeSize)) >= 1, name=f"sample_has_tumor[{s}]")
        else:
            # No surviving cluster has presence in this sample; allow root-only assignment
            m.addConstr(quicksum(omega[v,s] for v in range(1,treeSize)) == 0, name=f"sample_root_only[{s}]")

    # x,y
    x={}
    for v in range(treeSize):
        for s in sampleIDs:
            x[v,s]=m.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name=f"x[{v},{s}]")
    for s in sampleIDs:
        m.addConstr(quicksum(x[v,s] for v in range(treeSize))==1.0)

    # B: parent→child (edge) matrix only (A removed) — create only for allowed edges; fix others to 0
    B = {}
    for v in range(treeSize):
        for u in range(treeSize):
            if v == u or not allowed_edge[v][u]:
                B[v,u] = m.addVar(vtype=GRB.BINARY, lb=0, ub=0, name=f"B[{v},{u}]")  # fixed 0
            else:
                B[v,u] = m.addVar(vtype=GRB.BINARY, name=f"B[{v},{u}]")
    if fixed_B_edges is not None:
        for v in range(treeSize):
            for u in range(treeSize):
                if v == u:
                    continue
                m.addConstr(B[v, u] == (1 if (v, u) in fixed_B_edges else 0), name=f"fixed_B[{v},{u}]")

    # Enforce required ancestries from A_input using only B (path existence via transitive cover on A_input)
    if args.ancestryMatrixFile:
        for v in range(treeSize):
            for u in range(treeSize):
                if v != u and A_input[v][u] == 1:
                    m.addConstr(
                        (B[v,u] if (v,u) in B else 0) + quicksum((B[v,w] if (v,w) in B else 0) * A_input[w][u]
                                                                  for w in range(treeSize) if w != v) >= 1
                    )

    # tree structure on B
    for u in range(treeSize):
        if u == root:
            m.addConstr(quicksum(B[v,u] for v in range(treeSize) if v != u) == 0)
        else:
            m.addConstr(quicksum(B[v,u] for v in range(treeSize) if v != u) == 1)
    if treeSize > 1 and max_root_out_degree_arg is not None:
        m.addConstr(
            quicksum(B[root,u] for u in range(treeSize) if u != root) <= max_root_out_degree_arg,
            name="max_root_out_degree"
        )

    # y + linearization only across allowed edges
    y, q = {}, {}
    for v in range(treeSize):
        for s in sampleIDs:
            y[v,s] = m.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name=f"y[{v},{s}]")

    for v, u in ALLOWED_PAIRS:
        for s in sampleIDs:
            q[v,u,s] = m.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name=f"q[{v},{u},{s}]")
            m.addConstr(q[v,u,s] <= y[u,s])
            m.addConstr(q[v,u,s] <= B[v,u])
            m.addConstr(q[v,u,s] >= y[u,s] - (1 - B[v,u]))

    # Optional observed-child linearization for refined node trigger mode "proxy_all_obs":
    # q_obs[v,u,s] = B[v,u] * (sum_c CP[c,s] * delta[c,u]).
    q_obs = None
    if migration_presence_source == 'node' and node_trigger_mode == 'proxy_all_obs':
        q_obs = {}
        for v, u in ALLOWED_PAIRS:
            for s in sampleIDs:
                q_obs[v,u,s] = m.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=OBS_DCF_UB, name=f"q_obs[{v},{u},{s}]")
                obs_u_s = quicksum(CP[cid][s] * delta[cid, u] for cid in clusterIDs)
                m.addConstr(q_obs[v,u,s] <= obs_u_s)
                m.addConstr(q_obs[v,u,s] <= OBS_DCF_UB * B[v,u])
                m.addConstr(q_obs[v,u,s] >= obs_u_s - OBS_DCF_UB * (1 - B[v,u]))

    for v in range(treeSize):
        for s in sampleIDs:
            m.addConstr(y[v,s] == x[v,s] + quicksum(q[v2,u,s] for (v2,u) in ALLOWED_PAIRS if v2 == v))

    # Acyclicity via indicator constraints on depths (stronger than big-M)
    depth = {v: m.addVar(vtype=GRB.INTEGER, lb=0, ub=treeSize-1, name=f"depth[{v}]") for v in range(treeSize)}
    m.addConstr(depth[root] == 0)
    for v, u in ALLOWED_PAIRS:
        m.addGenConstrIndicator(B[v,u], True, depth[u] >= depth[v] + 1, name=f"acyc[{v},{u}]")

    # Anti-parallel cut to tighten LP: no 2-cycles
    for v in range(treeSize):
        for u in range(v+1, treeSize):
            m.addConstr(B[v,u] + B[u,v] <= 1)

    # Connectivity: single-commodity flow from root to all nodes (prevents forests and enforces reachability)
    BIG = treeSize - 1
    flow = {}
    for v, u in ALLOWED_PAIRS:
        flow[v,u] = m.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=BIG, name=f"flow[{v},{u}]")
        m.addConstr(flow[v,u] <= BIG * B[v,u])

    for u in range(treeSize):
        inflow  = quicksum(flow[v,u] for v in range(treeSize) if (v,u) in flow)
        outflow = quicksum(flow[u,w] for w in range(treeSize) if (u,w) in flow)
        if u == root:
            m.addConstr(outflow - inflow == treeSize - 1)
        else:
            m.addConstr(inflow - outflow == 1)

    # Redundant but tightening cut: total number of edges is N-1
    m.addConstr(quicksum(B[v,u] for v in range(treeSize) for u in range(treeSize) if v!=u) == treeSize - 1)

    # reward (optional) — only on allowed edges to keep model compact
    eq_presence={}
    match_all={}
    reward_expr=None
    duo_present_expr=None
    duo_absent_expr=None
    duo_presence=None
    if abs(co_presence_reward)>epsilon and ALLOWED_PAIRS:
        reward_expr = LinExpr()
        for (v,u) in ALLOWED_PAIRS:
            eqs=[]
            for s in sampleIDs:
                eq_presence[v,u,s]=m.addVar(vtype=GRB.BINARY, name=f"eq_presence[{v},{u},{s}]")
                # eq_presence[v,u,s] == 1  iff  omega[v,s] == omega[u,s]
                m.addConstr(eq_presence[v,u,s] <= 1 - omega[v,s] + omega[u,s])
                m.addConstr(eq_presence[v,u,s] <= 1 + omega[v,s] - omega[u,s])
                m.addConstr(eq_presence[v,u,s] >= omega[v,s] + omega[u,s] - 1)
                m.addConstr(eq_presence[v,u,s] >= 1 - omega[v,s] - omega[u,s])
                eqs.append(eq_presence[v,u,s])
            match_all[v,u]=m.addVar(vtype=GRB.BINARY, name=f"match_all[{v},{u}]")
            m.addConstr(match_all[v,u] <= B[v,u])
            for s in sampleIDs:
                m.addConstr(match_all[v,u] <= eq_presence[v,u,s])
            if sampleIDs:
                m.addConstr(match_all[v,u] >= B[v,u] + quicksum(eqs) - len(sampleIDs))
            reward_expr += match_all[v,u]
    if ALLOWED_PAIRS and (abs(duo_present_reward)>epsilon or abs(duo_absent_reward)>epsilon) and duo_presence_source != 'stacking':
        # Choose presence signal for duo rewards
        if duo_presence_source == 'observed':
            duo_presence = omega
        else:
            duo_presence = {}
            src = y if duo_presence_source == 'inferred' else x
            thr = duo_presence_threshold
            for v in range(treeSize):
                for s in sampleIDs:
                    duo_presence[v,s] = m.addVar(vtype=GRB.BINARY, name=f"duo_presence[{v},{s}]")
                    if v == root:
                        m.addConstr(duo_presence[v,s] == 0)
                        continue
                    # duo_presence == 1 iff src >= thr (within epsilon)
                    m.addConstr(src[v,s] >= thr * duo_presence[v,s])
                    m.addConstr(src[v,s] <= (thr - epsilon) + (1 - thr + epsilon) * duo_presence[v,s])
        if abs(duo_present_reward)>epsilon:
            duo_present_expr = LinExpr()
        if abs(duo_absent_reward)>epsilon:
            duo_absent_expr = LinExpr()
        for (v,u) in ALLOWED_PAIRS:
            for s in sampleIDs:
                if duo_present_expr is not None:
                    duo_present = m.addVar(vtype=GRB.BINARY, name=f"duo_present[{v},{u},{s}]")
                    m.addConstr(duo_present <= duo_presence[v,s])
                    m.addConstr(duo_present <= duo_presence[u,s])
                    m.addConstr(duo_present <= B[v,u])
                    m.addConstr(duo_present >= duo_presence[v,s] + duo_presence[u,s] + B[v,u] - 2)
                    duo_present_expr += duo_present
                if duo_absent_expr is not None:
                    duo_absent = m.addVar(vtype=GRB.BINARY, name=f"duo_absent[{v},{u},{s}]")
                    m.addConstr(duo_absent <= 1 - duo_presence[v,s])
                    m.addConstr(duo_absent <= 1 - duo_presence[u,s])
                    m.addConstr(duo_absent <= B[v,u])
                    m.addConstr(duo_absent >= B[v,u] - duo_presence[v,s] - duo_presence[u,s])
                    duo_absent_expr += duo_absent

    # Joint migration variables (main Section 2.3):
    # mu[v,site] indicates (v,site) is used in refined tree;
    # d[child,parent,site] indicates descent edge parent->child at site.
    mu_vars = None
    root_descent_vars = None
    descent_edge_vars = None
    omega_site_vars = None
    alpha_vars = None
    clone_site_pair_count_expr = None
    migration_count_expr = None
    if migration_enabled and site_list:
        mu_vars = {}
        root_descent_vars = None
        descent_edge_vars = {}
        def add_threshold_hit(signal_expr, hit_var, signal_lb, signal_ub, name):
            if migration_presence_threshold <= signal_lb + epsilon:
                m.addConstr(hit_var == 1, name=f"{name}_always")
            elif migration_presence_threshold > signal_ub + epsilon:
                m.addConstr(hit_var == 0, name=f"{name}_never")
            else:
                m.addConstr(
                    signal_expr >= signal_lb + (migration_presence_threshold - signal_lb) * hit_var,
                    name=f"{name}_lb"
                )
                m.addConstr(
                    signal_expr <= (migration_presence_threshold - epsilon) +
                                   (signal_ub - migration_presence_threshold + epsilon) * hit_var,
                    name=f"{name}_ub"
                )

        def sample_presence_signal(v, samp):
            signal_lb = 0.0
            signal_ub = 1.0
            if migration_presence_source == 'node':
                if node_trigger_mode in ('any_x', 'avg_x'):
                    signal_expr = x[v, samp]
                else:
                    node_obs = quicksum(CP[cid][samp] * delta[cid, v] for cid in clusterIDs)
                    if node_trigger_mode == 'proxy_all_obs':
                        child_mass = quicksum(q_obs[v2, u, samp] for (v2, u) in ALLOWED_PAIRS if v2 == v)
                    else:
                        child_mass = quicksum(q[v2, u, samp] for (v2, u) in ALLOWED_PAIRS if v2 == v)
                    signal_expr = node_obs - child_mass
                    signal_lb = -float(treeSize) * OBS_DCF_UB
                    signal_ub = OBS_DCF_UB
            elif migration_presence_source == 'lineage':
                signal_expr = y[v, samp]
            else:
                signal_expr = omega[v, samp]
            return signal_expr, signal_lb, signal_ub

        def add_site_presence_or(site_presence_var, v, site, samples_here, prefix, exact=True):
            hits = []
            for si, samp in enumerate(samples_here):
                hit = m.addVar(vtype=GRB.BINARY, name=f"{prefix}_sample_hit[{v},{site},{si}]")
                signal_expr, signal_lb, signal_ub = sample_presence_signal(v, samp)
                add_threshold_hit(signal_expr, hit, signal_lb, signal_ub, f"{prefix}_sample_hit[{v},{site},{si}]")
                m.addConstr(site_presence_var >= hit, name=f"{prefix}_site_or_lb[{v},{site},{si}]")
                hits.append(hit)
            if hits and exact:
                m.addConstr(site_presence_var <= quicksum(hits), name=f"{prefix}_site_or_ub[{v},{site}]")
            elif not hits:
                m.addConstr(site_presence_var == 0, name=f"{prefix}_site_no_samples[{v},{site}]")

        omega_site_vars = {}
        alpha_vars = {}
        alpha_child_product_vars = {}
        admissibility_bigM = treeSize

        for v in range(1, treeSize):
            for site in site_list:
                mu_vars[v, site] = m.addVar(vtype=GRB.BINARY, name=f"mu[{v},{site}]")
                omega_site_vars[v, site] = m.addVar(vtype=GRB.BINARY, name=f"omega_site[{v},{site}]")
                alpha_vars[v, site] = m.addVar(vtype=GRB.BINARY, name=f"alpha[{v},{site}]")

        for v, u in ALLOWED_PAIRS:
            if v == root or u == root:
                continue
            for site in site_list:
                alpha_child_product_vars[v, u, site] = m.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0,
                                                        name=f"alpha_child_product[{v},{u},{site}]")

        # For each non-root node, pick exactly one (parent,site) descent indicator.
        for child in range(treeSize):
            if child == root:
                continue
            for parent in range(treeSize):
                for site in site_list:
                    dvar = m.addVar(vtype=GRB.BINARY, name=f"origin_indicator_edge[{child},{parent},{site}]")
                    descent_edge_vars[child, parent, site] = dvar
                    m.addConstr(dvar <= mu_vars[child, site], name=f"migration_d_le_child[{child},{parent},{site}]")
                    if parent != root:
                        m.addConstr(dvar <= mu_vars[parent, site], name=f"migration_d_le_parent[{child},{parent},{site}]")
                    elif site in metastatic_site_set:
                        m.addConstr(dvar == 0, name=f"migration_root_child_not_metastatic_origin[{child},{site}]")
                    m.addConstr(dvar <= B[parent, child], name=f"migration_d_le_B[{child},{parent},{site}]")
            m.addConstr(
                quicksum(descent_edge_vars[child, parent, site] for parent in range(treeSize) for site in site_list) == 1,
                name=f"migration_d_onehot[{child}]"
            )

        # Paper admissibility:
        # - omega_site[v,site] is the observed local evidence (paper omega);
        # - alpha[v,site] is the admissibility indicator that recursively closes local
        #   evidence up the selected tree;
        # - mu[v,site] is then restricted to the admissible set via mu <= alpha.
        #
        # We also keep the "observed pairs must appear in R" behavior by forcing
        # mu >= omega_site. The healthy root has no site-specific omega/alpha/mu variables.
        for v in range(1, treeSize):
            for site in site_list:
                samples_here = site_to_samples.get(site, [])
                if not samples_here:
                    m.addConstr(mu_vars[v, site] == 0, name=f"migration_no_samples[{v},{site}]")
                    m.addConstr(omega_site_vars[v, site] == 0, name=f"omega_site_no_samples[{v},{site}]")
                    m.addConstr(alpha_vars[v, site] == 0, name=f"alpha_no_samples[{v},{site}]")
                    continue

                add_site_presence_or(omega_site_vars[v, site], v, site, samples_here, "omega_site")
                m.addConstr(mu_vars[v, site] >= omega_site_vars[v, site],
                            name=f"migration_obs_subset_m[{v},{site}]")

        for v, u in ALLOWED_PAIRS:
            if v == root or u == root:
                continue
            for site in site_list:
                m.addConstr(alpha_child_product_vars[v, u, site] <= B[v, u],
                            name=f"alpha_child_product_le_B[{v},{u},{site}]")
                m.addConstr(alpha_child_product_vars[v, u, site] <= alpha_vars[u, site],
                            name=f"alpha_child_product_le_childA[{v},{u},{site}]")
                m.addConstr(
                    alpha_child_product_vars[v, u, site] >= alpha_vars[u, site] + B[v, u] - 1,
                    name=f"alpha_child_product_ge_childA[{v},{u},{site}]"
                )

        for v in range(1, treeSize):
            outdeg_expr = quicksum(B[v, u] for u in range(treeSize) if (v, u) in B and v != u)
            for site in site_list:
                child_terms = quicksum(
                    alpha_child_product_vars[v, u, site]
                    for (vp, u) in ALLOWED_PAIRS
                    if vp == v and (v, u, site) in alpha_child_product_vars
                )
                m.addConstr(
                    alpha_vars[v, site] <= omega_site_vars[v, site] + child_terms,
                    name=f"alpha_recur_ub[{v},{site}]"
                )
                m.addConstr(
                    float(admissibility_bigM) * alpha_vars[v, site] >= omega_site_vars[v, site] + child_terms,
                    name=f"alpha_recur_lb[{v},{site}]"
                )
                m.addConstr(
                    mu_vars[v, site] <= alpha_vars[v, site],
                    name=f"migration_admissible_cap[{v},{site}]"
                )
            m.addConstr(
                quicksum(alpha_vars[v, site] for site in site_list) >= 1 - outdeg_expr,
                name=f"migration_leaf_admissible[{v}]"
            )

        clone_site_pair_count_expr = quicksum(mu_vars[v, site] for v in range(1, treeSize) for site in site_list)
        migration_count_expr = migration_count_expr_from_mu_vars(mu_vars, site_list, root_node=root)

        # Optional explicit count controls:
        # sum(mu[v,site]) <= M_pairs and/or counted migrations == M_mig.
        if max_clone_site_pairs_arg is not None:
            m.addConstr(
                clone_site_pair_count_expr <= float(max_clone_site_pairs_arg),
                name="num_clone_site_pairs_cap"
            )
        if num_migrations_arg is not None:
            m.addConstr(
                migration_count_expr == float(num_migrations_arg),
                name="migration_exact_migration_count"
            )

    if ALLOWED_PAIRS and duo_presence_source == 'stacking' and abs(duo_present_reward) > epsilon:
        if mu_vars is None:
            raise ValueError("--duo-presence-source stacking requires CITUPm migration site calls, but they were not built.")
        # PDF Section 7, generalized to site-level presence: \\\ depricated pdf location
        # p[v,u,site] = mu[v,site] AND mu[u,site] AND B[v,u],
        # with nonnegative edge surplus z[v,u] = sum_site p[v,u,site] - B[v,u].
        duo_present_expr = LinExpr()
        for (v, u) in ALLOWED_PAIRS:
            if v == root or u == root:
                continue
            duo_stack_sum = LinExpr()
            for site in site_list:
                duo_stack = m.addVar(vtype=GRB.BINARY, name=f"duo_stacking[{v},{u},{site}]")
                m.addConstr(duo_stack <= mu_vars[v, site], name=f"duo_stacking_le_tail[{v},{u},{site}]")
                m.addConstr(duo_stack <= mu_vars[u, site], name=f"duo_stacking_le_head[{v},{u},{site}]")
                m.addConstr(duo_stack <= B[v, u], name=f"duo_stacking_le_B[{v},{u},{site}]")
                m.addConstr(
                    duo_stack >= mu_vars[v, site] + mu_vars[u, site] + B[v, u] - 2,
                    name=f"duo_stacking_ge_and[{v},{u},{site}]"
                )
                duo_stack_sum += duo_stack
            duo_stack_surplus = m.addVar(vtype=GRB.INTEGER, lb=0.0, name=f"duo_stacking_surplus[{v},{u}]")
            m.addConstr(
                duo_stack_surplus == duo_stack_sum - B[v, u],
                name=f"duo_stacking_surplus_def[{v},{u}]"
            )
            duo_present_expr += duo_stack_surplus

    # fit objective
    residual_expr=None
    # Legacy bug mode: use the old mutation residual for BOTH objectives, so they match exactly
    if args.legacy_fit:
        residual_expr=QuadExpr()
        z={}
        def _mut_w(cid):
            k=float(clusterID_counts[cid])
            if args.clusterWeight=='equal': return 1.0/k if k>0 else 1.0
            if args.clusterWeight=='sqrt':  return (k**-0.5) if k>0 else 1.0
            return 1.0
        for i,mid in enumerate(mutIDs):
            cid=c[mid]
            for v in range(treeSize):
                for s in sampleIDs:
                    z[mid,v,s]=m.addVar(vtype=GRB.CONTINUOUS, lb=0.0, name=f"z[{mid},{v},{s}]")
                    # Legacy fit: implicit M=1 on non-selected nodes
                    m.addConstr(z[mid,v,s] >= (delta[cid,v] - 1) + (f[mid][s] - y[v,s]))
                    m.addConstr(z[mid,v,s] >= (delta[cid,v] - 1) - (f[mid][s] - y[v,s]))
                    residual_expr += _mut_w(cid) * z[mid,v,s] * z[mid,v,s]
    elif args.objective=='mutation':
        residual_expr=QuadExpr()
        z={}
        def _mut_w(cid):
            k=float(clusterID_counts[cid])
            if args.clusterWeight=='equal': return 1.0/k if k>0 else 1.0
            if args.clusterWeight=='sqrt':  return (k**-0.5) if k>0 else 1.0
            return 1.0
        for i,mid in enumerate(mutIDs):
            cid=c[mid]
            if args.fixClusterMapping:
                # Mapping is fixed: only the fixed node can host this cluster.
                vfix = clusterIndex[cid]
                for s in sampleIDs:
                    z[mid,vfix,s]=m.addVar(vtype=GRB.CONTINUOUS, lb=0.0, name=f"z[{mid},{vfix},{s}]")
                    # Exact absolute value without big-M, since delta[cid,vfix]==1
                    m.addConstr(z[mid,vfix,s] >=  (f[mid][s] - y[vfix,s]))
                    m.addConstr(z[mid,vfix,s] >= -(f[mid][s] - y[vfix,s]))
                    residual_expr += _mut_w(cid) * z[mid,vfix,s] * z[mid,vfix,s]
            else:
                for v in range(treeSize):
                    for s in sampleIDs:
                        z[mid,v,s]=m.addVar(vtype=GRB.CONTINUOUS, lb=0.0, name=f"z[{mid},{v},{s}]")
                        # Use a safe big-M so when delta==0 these bounds are non-binding (<=0)
                        m.addConstr(z[mid,v,s] >= BIGM_MUT_ABS*(delta[cid,v] - 1) + (f[mid][s] - y[v,s]))
                        m.addConstr(z[mid,v,s] >= BIGM_MUT_ABS*(delta[cid,v] - 1) - (f[mid][s] - y[v,s]))
                        residual_expr += _mut_w(cid) * z[mid,v,s] * z[mid,v,s]
    else:
        w,t={},{}
        for cid in clusterIDs:
            for v in range(treeSize):
                for s in sampleIDs:
                    w[cid,v,s]=m.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name=f"w[{cid},{v},{s}]")
                    m.addConstr(w[cid,v,s] <= y[v,s])
                    m.addConstr(w[cid,v,s] <= delta[cid,v])
                    m.addConstr(w[cid,v,s] >= y[v,s] - (1 - delta[cid,v]))
            for s in sampleIDs:
                t[cid,s]=m.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name=f"t[{cid},{s}]")
                m.addConstr(t[cid,s] == quicksum(w[cid,v,s] for v in range(treeSize)))
        residual_expr=QuadExpr()
        for cid in clusterIDs:
            k=float(clusterID_counts[cid])
            wgt = 1.0 if args.clusterWeight=='equal' else (k**0.5 if args.clusterWeight=='sqrt' else k)
            for s in sampleIDs:
                residual_expr += wgt * (t[cid,s]-CP[cid][s])*(t[cid,s]-CP[cid][s])

    fit_obj = residual_expr + 0
    if reward_expr is not None:
        fit_obj -= co_presence_reward * reward_expr
    if duo_present_expr is not None:
        fit_obj -= duo_present_reward * duo_present_expr
    if duo_absent_expr is not None:
        fit_obj -= duo_absent_reward * duo_absent_expr
    if clone_site_pair_count_expr is not None and abs(migration_weight) > epsilon:
        fit_obj += migration_weight * migration_count_expr

    bounded_score_expr = residual_expr if args.migration_score_source == 'without-reward' else fit_obj
    if score_bounded_migrations:
        m.addQConstr(bounded_score_expr <= score_upper, name='migration_score_cap')

    return dict(model=m, root=root, delta=delta, omega=omega, x=x, y=y, B=B,
                fit_obj=fit_obj, residual_expr=residual_expr, reward_expr=reward_expr,
                duo_present_expr=duo_present_expr, duo_absent_expr=duo_absent_expr,
                mu_vars=mu_vars, root_descent_vars=root_descent_vars,
                descent_edge_vars=descent_edge_vars, omega_site_vars=omega_site_vars,
                alpha_vars=alpha_vars, clone_site_pair_count_expr=clone_site_pair_count_expr,
                migration_count_expr=migration_count_expr,
                bounded_score_expr=bounded_score_expr)

# ----------------- Solve pipeline -----------------
# storage
solutions = []  # collected trees: {core, b_edges, fit, residual_stage1, reward_stage1, stage1_delta, stage1_x, stage1_y}
enumerated_b = []  # all tree structures already enumerated (feasible or not)

best_fit = None
fit_upper = None
migration_minimum = None
migration_minimum_proven = False
score_tiebreak_proven = False

# The optional mode bounds the score and minimizes migration count first.
core = build_core()
m = core["model"]
if score_bounded_migrations:
    m.setObjective(core['migration_count_expr'], GRB.MINIMIZE)
    m.Params.MIPGap = 0
    m.Params.MIPGapAbs = 0
    print(f'[info] Minimizing migrations with {args.migration_score_source} score '
          f'<= {score_upper:.17g} (reference {score_reference:.17g}).')
    apply_first_phase_warm_start(core, warm_start_B_edges, warm_start_delta_assign)
else:
    m.setObjective(core["fit_obj"], GRB.MINIMIZE)
    m.Params.MIPGap    = float(args.fit_gap)
    m.Params.MIPGapAbs = float(args.fit_gap_abs)
if args.fit_time_limit is not None: m.Params.TimeLimit = float(args.fit_time_limit)
if args.write_mps:
    m.update()
    m.write(args.write_mps)
    print(f"[info] Wrote stage-1 model to {args.write_mps}")
if args.write_mps_only:
    if not args.write_mps:
        raise ValueError("--write-mps-only requires --write-mps PATH")
    print("[info] --write-mps-only requested; exiting before stage-1 optimization.")
    sys.exit(0)
stage1_started = time.monotonic()
m.optimize()
if m.Status in (GRB.INFEASIBLE, GRB.UNBOUNDED, GRB.INF_OR_UNBD):
    if getattr(args, 'debug_iis', False):
        try:
            m.computeIIS()
            m.write(args.pathOutputFilePrefix + ".stage1_best.iis.ilp")
            m.write(args.pathOutputFilePrefix + ".stage1_best.lp")
        except: pass
    raise RuntimeError("Stage-1 infeasible/unbounded.")
if m.SolCount==0: raise RuntimeError("Stage-1: no incumbent.")

# Enforce optimality unless explicitly allowed (skip for resource limits)
allow_stage1_subopt = (m.Status in (GRB.TIME_LIMIT, GRB.MEM_LIMIT))
require_opt_stage1 = getattr(args, 'require_optimal_stage1', True)
if m.Status != GRB.OPTIMAL and require_opt_stage1 and not allow_stage1_subopt:
    try:
        bound = m.ObjBound
        val = m.ObjVal
        gap = m.MIPGap
        print(f"[error] Stage-1 not proven optimal (status={m.Status}). Incumbent={val:.6g}, bestBound={bound:.6g}, gap={gap:.3g}", file=sys.stderr)
    except Exception:
        print(f"[error] Stage-1 not proven optimal (status={m.Status}).", file=sys.stderr)
    raise RuntimeError("Stage-1 did not finish optimal and --require-optimal-stage1 is set. Use --allow-suboptimal-stage1 to proceed anyway.")
elif m.Status != GRB.OPTIMAL and allow_stage1_subopt:
    try:
        val = m.ObjVal
        gap = m.MIPGap
        print(f"[warn] Stage-1 ended with status={m.Status}; continuing with incumbent. Incumbent={val:.6g}, gap={gap:.3g}", file=sys.stderr)
    except Exception:
        print(f"[warn] Stage-1 ended with status={m.Status}; continuing with incumbent.", file=sys.stderr)

if score_bounded_migrations:
    migration_primary_status = m.Status
    migration_primary_bound = m.ObjBound
    migration_minimum = roundInt(core['migration_count_expr'].getValue())
    migration_minimum_proven = m.Status == GRB.OPTIMAL
    if not migration_minimum_proven:
        print('[warn] Minimum migration count is not proven; reporting the best '
              'incumbent count within the score bound.', file=sys.stderr)

    if args.migration_score_tiebreak:
        # Fix the primary objective before choosing the best score among its ties.
        for var in m.getVars():
            var.Start = var.X
        m.addConstr(core['migration_count_expr'] == migration_minimum,
                    name='minimum_migration_count')
        stage1_limit = args.fit_time_limit if args.fit_time_limit is not None else args.gurobi_max_time
        remaining = None if stage1_limit is None else float(stage1_limit) - (time.monotonic() - stage1_started)
        score_tiebreak_status = 'not_run_budget_exhausted'
        if remaining is None or remaining > 0:
            m.setObjective(core['bounded_score_expr'], GRB.MINIMIZE)
            if remaining is not None:
                m.Params.TimeLimit = remaining
            m.optimize()
            if m.SolCount == 0:
                raise RuntimeError('Score tie-breaking found no incumbent at the feasible '
                                   f'migration count {migration_minimum} (status={m.Status}).')
            if m.Status not in (GRB.OPTIMAL, GRB.TIME_LIMIT, GRB.MEM_LIMIT):
                raise RuntimeError(f'Score tie-breaking failed (status={m.Status}).')
            score_tiebreak_status = str(m.Status)
            score_tiebreak_proven = m.Status == GRB.OPTIMAL
    else:
        score_tiebreak_status = 'not_run_disabled'
        print('[info] Internal score tie-break disabled; keeping the migration incumbent '
              'and continuing migration/comigration refinement.')

    attained_score = core['bounded_score_expr'].getValue()
    print(f'[info] Migration count {migration_minimum}; minimum proven '
          f'{migration_minimum_proven}; bounded score {attained_score:.17g}.')
    with open(args.pathOutputFilePrefix + '.migration_score_optimization.tsv', 'w') as report:
        report.write('key\tvalue\n')
        for key, value in (
            ('score_source', args.migration_score_source),
            ('score_reference', format(score_reference, '.17g')),
            ('tolerance_mode', 'additive' if args.migration_score_abs_tolerance is not None else 'paper_eq9'),
            ('tolerance', format(tolerance, '.17g')),
            ('score_upper_bound', format(score_upper, '.17g')),
            ('migration_count', migration_minimum),
            ('migration_minimum_proven', 'YES' if migration_minimum_proven else 'NO'),
            ('migration_solver_status', migration_primary_status),
            ('migration_objective_bound', format(migration_primary_bound, '.17g')),
            ('score_tiebreak_enabled', 'YES' if args.migration_score_tiebreak else 'NO'),
            ('score_tiebreak_status', score_tiebreak_status),
            ('score_tiebreak_proven', 'YES' if score_tiebreak_proven else 'NO'),
            ('attained_score', format(attained_score, '.17g')),
        ):
            report.write(f'{key}\t{value}\n')

best_fit = core["fit_obj"].getValue()
# Base enumeration settings (near-optimal by default)
fit_upper = score_upper if score_bounded_migrations else best_fit + max(float(args.enum_fit_abs), float(args.enum_fit_rel)*abs(best_fit))
all_pairs = set(ALLOWED_PAIRS)
enum_deadline = time.time() + float(args.enum_time_limit) if args.enum_time_limit is not None else None
_max_t = int(args.enum_max_trees)
_unlimited = (_max_t <= 0)
target_solutions = float("inf") if _unlimited else _max_t
# K-best mode: enumerate the K best-fit trees (in increasing fit), regardless of ties
kbest_mode = False
if getattr(args, 'num_solutions', None) not in (None, 0):
    _max_t = max(1, int(args.num_solutions))
    _unlimited = False
    target_solutions = _max_t
    kbest_mode = True  # do not cap fit by best_fit; let next solves improve objective to next-best

def capture_B(core):
    B = core['B']
    edges = set()
    for v in range(treeSize):
        for u in range(treeSize):
            if v != u and roundVar(B[v, u]) == 1:
                edges.add((v, u))
    return edges


def capture_stage1_metrics(core):
    fit_val = core["fit_obj"].getValue()
    residual_val = core["residual_expr"].getValue() if core["residual_expr"] is not None else fit_val
    reward_val = 0.0
    if core.get("reward_expr") is not None:
        reward_val += co_presence_reward * core["reward_expr"].getValue()
    if core.get("duo_present_expr") is not None:
        reward_val += duo_present_reward * core["duo_present_expr"].getValue()
    if core.get("duo_absent_expr") is not None:
        reward_val += duo_absent_reward * core["duo_absent_expr"].getValue()
    return fit_val, residual_val, reward_val

def snapshot_stage1_solution(core):
    delta_assign = {}
    for cid in clusterIDs:
        assigned = None
        for v in range(treeSize):
            if roundInt(core["delta"][cid,v].X) == 1:
                assigned = v
                break
        delta_assign[cid] = assigned
    x_vals = {(v,s): core["x"][v,s].X for v in range(treeSize) for s in sampleIDs}
    y_vals = {(v,s): core["y"][v,s].X for v in range(treeSize) for s in sampleIDs}
    return delta_assign, x_vals, y_vals

def compute_root_degree(sol):
    edges = sol.get("b_edges")
    if not edges:
        edges = set()
        B = sol["core"]["B"]
        for v in range(treeSize):
            for u in range(treeSize):
                if v != u and roundVar(B[v, u]) == 1:
                    edges.add((v, u))
    return sum(1 for (p, ch) in edges if p == 0)

def write_solution(sol_idx, sol, tag=""):
    core_i = sol["core"]
    sidecar_suffix = f"_sol{sol_idx}{tag}"
    fit_val = sol.get("fit")
    if fit_val is None:
        try:
            fit_val = core_i["fit_obj"].getValue()
        except Exception:
            fit_val = 0.0
    residual_val = sol.get("residual_stage1")
    if residual_val is None:
        try:
            residual_val = core_i["residual_expr"].getValue() if core_i["residual_expr"] is not None else fit_val
        except Exception:
            residual_val = fit_val
    reward_val = sol.get("reward_stage1")
    if reward_val is None:
        try:
            reward_val = 0.0
            if core_i.get("reward_expr") is not None:
                reward_val += co_presence_reward * core_i["reward_expr"].getValue()
            if core_i.get("duo_present_expr") is not None:
                reward_val += duo_present_reward * core_i["duo_present_expr"].getValue()
            if core_i.get("duo_absent_expr") is not None:
                reward_val += duo_absent_reward * core_i["duo_absent_expr"].getValue()
        except Exception:
            reward_val = 0.0
    migration_status = sol.get("migration_status")
    b_edges = sorted(sol.get("b_edges", []))

    stage1_delta = sol.get("stage1_delta") or {}
    stage1_x = sol.get("stage1_x") or {}
    stage1_y = sol.get("stage1_y") or {}

    # cluster↦node
    nodeOfCluster = {}
    for cid in clusterIDs:
        node = None
        for v in range(treeSize):
            if roundVar(core_i["delta"][cid, v]) == 1:
                node = v
                break
        nodeOfCluster[cid] = node

    edges = []
    for v in range(treeSize):
        for u in range(treeSize):
            if v != u and roundVar(core_i["B"][v, u]) == 1:
                edges.append((v, u))
    if not edges:
        edges = b_edges

    rootDegree = sum(1 for (p, ch) in edges if p == 0)

    # A migration-only solve does not certify minimum fit at the chosen count.
    fit_proven = (core_i["model"].Status == GRB.OPTIMAL and
                  (not score_bounded_migrations or sol_idx > 0 or score_tiebreak_proven))

    # tree file
    with open(args.pathOutputFilePrefix + f".tree_{sol_idx}{tag}.txt", "w") as fp:
        fp.write("TREE_SCORE\n" + format(fit_val, '.17g') + "\n")
        fp.write("TREE_SCORE_NO_REWARD\n" + format(residual_val, '.17g') + "\n")
        fp.write("REWARD_CONTRIBUTION\n" + format(reward_val, '.17g') + "\n")
        if score_bounded_migrations:
            fp.write('OPTIMIZATION_MODE\nminimize_migrations_at_score\n')
            fp.write('MIGRATION_MINIMUM_PROVEN\n' + ('YES' if migration_minimum_proven else 'NO') + '\n')
            fp.write('MIGRATION_SCORE_TIEBREAK_ENABLED\n' + ('YES' if args.migration_score_tiebreak else 'NO') + '\n')
            fp.write('MIGRATION_SCORE_TIEBREAK_STATUS\n' + score_tiebreak_status + '\n')
            fp.write('MIGRATION_SCORE_TIEBREAK_PROVEN\n' + ('YES' if score_tiebreak_proven else 'NO') + '\n')
        fp.write("NUM_MIGRATIONS\n" + ("NA" if sol.get("num_migrations") is None else str(sol["num_migrations"])) + "\n")
        fp.write("NUM_COMIGRATIONS\n" + ("NA" if sol.get("num_comigrations") is None else str(sol["num_comigrations"])) + "\n")
        if sol.get("refined_solution_index") is not None:
            fp.write("REFINED_SOLUTION_INDEX\n" + str(sol.get("refined_solution_index")) + "\n")
        if sol.get("refined_solution_count") is not None:
            fp.write("REFINED_SOLUTION_COUNT\n" + str(sol.get("refined_solution_count")) + "\n")
        if sol.get("refined_enumeration_status") is not None:
            fp.write("REFINED_ENUMERATION_STATUS\n" + str(sol.get("refined_enumeration_status")) + "\n")
        if sol.get("refined_enumeration_exhaustive") is not None:
            fp.write("REFINED_ENUMERATION_EXHAUSTIVE\n" + ("YES" if sol.get("refined_enumeration_exhaustive") else "NO") + "\n")
        fp.write("NUM_CLONE_SITE_PAIRS\n" + ("NA" if sol.get("num_clone_site_pairs") is None else str(sol["num_clone_site_pairs"])) + "\n")
        fp.write("HEALTHY_ROOT_NO_SITE\nYES\n")
        fp.write("IGNORE_ROOT_MIGRATIONS\nYES\n")
        fp.write("METASTATIC_SITES\n" + ("NA" if not metastatic_sites else ",".join(metastatic_sites)) + "\n")
        if migration_status:
            fp.write("MIGRATION_STATUS\n" + migration_status + "\n")
        if sol.get("migration_dag_before_postprocess") is not None:
            fp.write("SITE_DAG_BEFORE_POSTPROCESS\n" + str(sol.get("migration_dag_before_postprocess")) + "\n")
        if sol.get("migration_dag_postprocess_status") is not None:
            fp.write("DAG_POSTPROCESS_STATUS\n" + str(sol.get("migration_dag_postprocess_status")) + "\n")
        if sol.get("migration_labeling_source") is not None:
            fp.write("LABELING_SOURCE\n" + str(sol.get("migration_labeling_source")) + "\n")
        if sol.get("migration_relabel_change_count") is not None:
            fp.write("RELABEL_CHANGE_COUNT\n" + str(sol.get("migration_relabel_change_count")) + "\n")
        if sol.get("migration_dag_enum_status") is not None:
            fp.write("DAG_ENUM_STATUS\n" + str(sol.get("migration_dag_enum_status")) + "\n")
        if sol.get("migration_dag_enum_objective") is not None:
            fp.write("DAG_ENUM_OBJECTIVE\n" + str(sol.get("migration_dag_enum_objective")) + "\n")
        if sol.get("migration_dag_solutions") is not None:
            fp.write("DAG_COLORING_COUNT\n" + str(len(sol.get("migration_dag_solutions") or [])) + "\n")
        fp.write("SITE_DAG_FINAL\n" + migration_graph_dag_flag(sol) + "\n")
        fp.write("FIT_SOLVER_STATUS\n" + str(core_i["model"].Status) + "\n")
        fp.write("FIT_OPTIMAL\n" + ("YES" if fit_proven else "NO") + "\n")
        fp.write("ROOT_DEGREE\n" + str(rootDegree) + "\n")
        fp.write("\nPARENT_NODE\tCHILD_NODE\n")
        for (p, ch) in edges:
            fp.write(f"{p}\t{ch}\n")
        fp.write("\nCLUSTER\tNODE_OF_CLUSTER\n")
        for cid in clusterIDs:
            fp.write(f"{cid}\t{nodeOfCluster[cid]}\n")
        fp.write("\nMUTATION\tNODE_OF_MUTATION\n")
        for mid in mutIDs:
            fp.write(f"{mid}\t{nodeOfCluster.get(c[mid])}\n")
        fp.write("\nNODE\tNODE_FREQS\tCLADE_FREQS\n")
        for v in range(treeSize):
            node_vals = [stage1_x.get((v, s), safe_val(core_i["x"][v, s])) for s in sampleIDs]
            clade_vals = [stage1_y.get((v, s), safe_val(core_i["y"][v, s])) for s in sampleIDs]
            fp.write(str(v) + "\t" + ",".join(floatToStr(val + epsilon) for val in node_vals) + "\t")
            fp.write(",".join(floatToStr(val + epsilon) for val in clade_vals) + "\n")
        fp.write("\nMUT_ID\tOBSERVED_FREQS\tINFERRED_FREQS\tABS_OBSERVED_MINUS_INFERRED\n")
        for mid in mutIDs:
            v = nodeOfCluster.get(c[mid])
            fp.write(mid + "\t" + ",".join(floatToStr(f[mid][s]) for s in sampleIDs) + "\t")
            if v is None:
                inferred = [0.0 for _ in sampleIDs]
            else:
                inferred = [stage1_y.get((v, s), safe_val(core_i["y"][v, s])) for s in sampleIDs]
            fp.write(",".join(floatToStr(val) for val in inferred) + "\t")
            fp.write(",".join(floatToStr(abs(f[mid][s] - inferred[idx])) for idx, s in enumerate(sampleIDs)) + "\n")
        if sol.get("refined_nodes"):
            fp.write("\nREFINED_NODE\tSITE\tPRESENCE_SIGNAL\tOBSERVED\tMU\tD\tIS_ORIGIN\n")
            for rn in sol["refined_nodes"]:
                fp.write(f"{rn['node']}\t{rn['site']}\t{floatToStr(rn['presence_signal'])}\t{rn['observed']}\t{rn['mu']}\t{rn['d']}\t{rn['is_origin']}\n")
        if sol.get("mu"):
            fp.write("\nMU_ASSIGNMENT\tSITE\tMU\n")
            for v in range(1, treeSize):
                for site in site_list:
                    fp.write(f"{v}\t{site}\t{int(sol['mu'].get((v, site), 0))}\n")
        if sol.get("origin_site"):
            fp.write("\nORIGIN_SITE\nNODE\tSITE\n")
            fp.write("0\tNA_HEALTHY_ROOT_NO_SITE\n")
            for v in range(1, treeSize):
                fp.write(f"{v}\t{sol['origin_site'].get(v)}\n")
        if sol.get("migration_dag_solutions"):
            fp.write("\nDAG_COLORING\tSAME_ORIGIN_SCORE\tMIGRATION_EDGE_COUNT\tDAG\n")
            for dag_sol in sol["migration_dag_solutions"]:
                fp.write(f"{dag_sol.get('index')}\t{dag_sol.get('same_origin_score')}\t{dag_sol.get('migration_edge_count')}\t{dag_sol.get('dag')}\n")
            fp.write("\nDAG_COLORING_INDEX\tNODE\tSITE\n")
            for dag_sol in sol["migration_dag_solutions"]:
                for v in range(1, treeSize):
                    fp.write(f"{dag_sol.get('index')}\t{v}\t{dag_sol.get('origin', {}).get(v)}\n")
        if sol.get("migration_location_graph_edges") is not None:
            fp.write("\nLOCATION_FROM_SITE\tLOCATION_TO_SITE\n")
            for source, destination in sol.get("migration_location_graph_edges") or []:
                fp.write(f"{source}\t{destination}\n")
        if sol.get("refined_edges"):
            fp.write("\nREFINED_FROM_NODE\tREFINED_FROM_SITE\tREFINED_TO_NODE\tREFINED_TO_SITE\tEDGE_TYPE\n")
            for (fn, fs, tn, ts, etype) in sol["refined_edges"]:
                fp.write(f"{fn}\t{fs}\t{tn}\t{ts}\t{etype}\n")
            fp.write("\nMIGRATION_FROM_NODE\tMIGRATION_FROM_SITE\tMIGRATION_TO_NODE\tMIGRATION_TO_SITE\n")
            for (fn, fs, tn, ts, etype) in sol["refined_edges"]:
                if etype == "migration":
                    fp.write(f"{fn}\t{fs}\t{tn}\t{ts}\n")

    # Additional CITUPm migration outputs for easier downstream use
    if sol.get("refined_nodes"):
        with open(args.pathOutputFilePrefix + f".refined_nodes{sidecar_suffix}.tsv", "w") as fn:
            fn.write("node\tsite\tpresence_signal\tobserved\tmu\td\tis_origin\n")
            for rn in sol["refined_nodes"]:
                fn.write(f"{rn['node']}\t{rn['site']}\t{floatToStr(rn['presence_signal'])}\t{rn['observed']}\t{rn['mu']}\t{rn['d']}\t{rn['is_origin']}\n")
    if sol.get("refined_edges"):
        with open(args.pathOutputFilePrefix + f".refined_edges{sidecar_suffix}.tsv", "w") as fe:
            fe.write("from_node\tfrom_site\tto_node\tto_site\ttype\n")
            for (fn, fs, tn, ts, etype) in sol["refined_edges"]:
                fe.write(f"{fn}\t{fs}\t{tn}\t{ts}\t{etype}\n")
    if sol.get("migration_dag_solutions"):
        with open(args.pathOutputFilePrefix + f".migration_dag_colorings{sidecar_suffix}.tsv", "w") as fd:
            fd.write("coloring_index\tnode\tsite\tsame_origin_score\tmigration_edge_count\tdag\n")
            for dag_sol in sol["migration_dag_solutions"]:
                for v in range(1, treeSize):
                    fd.write(f"{dag_sol.get('index')}\t{v}\t{dag_sol.get('origin', {}).get(v)}\t{dag_sol.get('same_origin_score')}\t{dag_sol.get('migration_edge_count')}\t{dag_sol.get('dag')}\n")
    if sol.get("refined_solutions"):
        refined_solutions = sol.get("refined_solutions") or []
        with open(args.pathOutputFilePrefix + f".refined_solutions{sidecar_suffix}.tsv", "w") as fr:
            fr.write("refined_solution_index\tnode\tsite\tmu\tis_origin\n")
            for refined in refined_solutions:
                ridx = refined.get("refined_solution_index")
                m_map = refined.get("mu") or {}
                origin = refined.get("origin_site") or {}
                for v in range(1, treeSize):
                    for site in site_list:
                        is_origin = 1 if origin.get(v) == site else 0
                        fr.write(f"{ridx}\t{v}\t{site}\t{int(m_map.get((v, site), 0))}\t{is_origin}\n")
        with open(args.pathOutputFilePrefix + f".refined_solution_summary{sidecar_suffix}.tsv", "w") as frs:
            frs.write("refined_solution_index\tnum_clone_site_pairs\tnum_migrations\tnum_comigrations\tlocation_graph_edges\tDAG\n")
            for refined in refined_solutions:
                ridx = refined.get("refined_solution_index")
                arcs = ",".join(f"{source}->{destination}" for source, destination in refined.get("migration_location_graph_edges") or [])
                frs.write(f"{ridx}\t{refined.get('num_clone_site_pairs', 'NA')}\t{refined.get('num_migrations', 'NA')}\t{refined.get('num_comigrations', 'NA')}\t{arcs}\t{migration_graph_dag_flag(refined)}\n")
        with open(args.pathOutputFilePrefix + f".location_graphs{sidecar_suffix}.tsv", "w") as flg:
            flg.write("refined_solution_index\tfrom_site\tto_site\n")
            for refined in refined_solutions:
                ridx = refined.get("refined_solution_index")
                for source, destination in refined.get("migration_location_graph_edges") or []:
                    flg.write(f"{ridx}\t{source}\t{destination}\n")

def capture_integrated_migration(sol_idx, sol):
    core_i = sol["core"]
    root_node = core_i.get("root", 0)
    mu_vars = core_i.get("mu_vars")
    if (not migration_enabled) or mu_vars is None:
        return dict(migration_status="skipped", num_migrations=None,
                    num_clone_site_pairs=None, num_comigrations=None)

    edge_set = sol.get("b_edges", set())
    if not edge_set:
        edge_set = capture_B(core_i)
    parent_map, _ = build_tree_maps(edge_set)
    site_signal = compute_site_signal(sol, site_to_samples, parent_map=parent_map)
    _, observed = build_mu_assignment_bounds(site_signal, parent_map, site_list,
                                             root_node=root_node, core=core_i)

    incumbent_m = {(v, site): roundVar(mu_vars[v, site])
                   for v in range(1, treeSize) for site in site_list}
    target_num_migrations = migration_count_from_mu(
        incumbent_m, site_list, root_node=root_node
    )
    enum_result = enumerate_min_comigration_refined_trees(
        sol_idx, sol, parent_map, site_list, site_signal, observed,
        target_num_migrations, max_solutions=max_refined_solutions_arg
    )
    refined_solutions = enum_result.get("solutions") or []
    if not refined_solutions:
        return build_core_migration_fallback(
            sol, parent_map, site_list, site_signal, observed,
            target_num_migrations, enum_result, root_node=root_node
        )

    for idx, refined in enumerate(refined_solutions):
        refined["refined_solution_index"] = idx
        refined["refined_solution_count"] = len(refined_solutions)
        refined["refined_enumeration_status"] = enum_result.get("status", "unknown")
        refined["refined_enumeration_exhaustive"] = enum_result.get("exhaustive", False)

    primary = refined_solutions[0]
    primary.update(
        refined_solutions=refined_solutions
    )
    return primary

def finalize_solution(sol_idx, sol):
    sol["root_degree"] = compute_root_degree(sol)
    if migration_enabled:
        try:
            res = capture_integrated_migration(sol_idx, sol)
            sol.update(res)
            if res.get("migration_status") != "ok":
                print(f"[warn] Integrated migration extraction not ok for tree {sol_idx}: {res.get('migration_status')}", file=sys.stderr)
            if not res.get("refined_enumeration_exhaustive", True):
                print(f"[warn] Refined-tree enumeration for tree {sol_idx} was not exhaustive: {res.get('refined_enumeration_status')}", file=sys.stderr)
        except Exception as e:
            sol["migration_status"] = "error"
            sol["migration_error"] = str(e)
            sol["num_migrations"] = None
            sol["num_clone_site_pairs"] = None
            sol["num_comigrations"] = None
            print(f"[warn] Integrated migration extraction error for tree {sol_idx}: {e}", file=sys.stderr)
    else:
        sol["migration_status"] = "skipped"
    write_solution(sol_idx, sol, tag="")
    for alt in sol.get("refined_solutions", [])[1:]:
        alt_sol = dict(sol)
        alt_sol.update(alt)
        alt_sol.pop("refined_solutions", None)
        alt_idx = alt.get("refined_solution_index")
        write_solution(sol_idx, alt_sol, tag=f"_refined_{alt_idx}")

nogood_counter = 0
def add_hamming_cut(model, B_vars, edge_set, tag):
    global nogood_counter
    nogood_counter += 1
    model.addConstr(
        quicksum((1 - B_vars[v,u]) for (v,u) in edge_set) +
        quicksum(B_vars[v,u] for (v,u) in (all_pairs - edge_set))
        >= 1,
        name=f"{tag}_{nogood_counter}"
    )

def enumerate_next(prev_edge_sets, label_prefix):
    core_k = build_core()
    mk = core_k["model"]
    # In K-best mode, do NOT cap by fit_upper; otherwise, restrict to near-optimal window
    if score_bounded_migrations:
        mk.addConstr(core_k['migration_count_expr'] == migration_minimum,
                     name='minimum_migration_count')
    elif not kbest_mode:
        mk.addQConstr(core_k["fit_obj"] <= fit_upper, name="enum_fit_cap")
    for b_prev in prev_edge_sets:
        if not b_prev:
            continue
        add_hamming_cut(mk, core_k["B"], b_prev, label_prefix)
    enumeration_score_expr = core_k['bounded_score_expr'] if score_bounded_migrations else core_k['fit_obj']
    mk.setObjective(enumeration_score_expr, GRB.MINIMIZE)
    mk.Params.MIPGap = 0; mk.Params.MIPGapAbs = 0
    if args.fit_time_limit is not None:
        mk.Params.TimeLimit = float(args.fit_time_limit)

    attempts = 0
    while True:
        if enum_deadline is not None and time.time() >= enum_deadline:
            return None
        if score_bounded_migrations and enum_deadline is not None:
            mk.Params.TimeLimit = min(mk.Params.TimeLimit, enum_deadline - time.time())
        mk.optimize()
        if mk.Status not in (GRB.OPTIMAL, GRB.TIME_LIMIT) or mk.SolCount == 0:
            return None
        fit_k = core_k["fit_obj"].getValue()
        bounded_score_k = enumeration_score_expr.getValue()
        if score_bounded_migrations and bounded_score_k > score_upper + mk.Params.FeasibilityTol:
            return None
        if (not score_bounded_migrations) and (not kbest_mode) and (fit_k > fit_upper + 1e-12):
            return None
        b_edges_k = capture_B(core_k)
        if any(b_edges_k == b_prev for b_prev in prev_edge_sets):
            add_hamming_cut(mk, core_k["B"], b_edges_k, f"{label_prefix}_retry")
            mk.update()
            attempts += 1
            if attempts >= 5:
                return None
            continue
        f0, r0, rw0 = capture_stage1_metrics(core_k)
        delta_snap, x_snap, y_snap = snapshot_stage1_solution(core_k)
        return dict(core=core_k, b_edges=b_edges_k, fit=f0, residual_stage1=r0, reward_stage1=rw0, stage1_delta=delta_snap, stage1_x=x_snap, stage1_y=y_snap)

# Seed enumeration with best-fit tree
f0, r0, rw0 = capture_stage1_metrics(core)
delta_snap, x_snap, y_snap = snapshot_stage1_solution(core)
initial_sol = dict(core=core, b_edges=capture_B(core), fit=f0, residual_stage1=r0, reward_stage1=rw0, stage1_delta=delta_snap, stage1_x=x_snap, stage1_y=y_snap)
enumerated_b.append(initial_sol["b_edges"])
solutions.append(initial_sol)
finalize_solution(0, initial_sol)

# Enumerate near-optimal trees until target reached or infeasible
while (_unlimited or len(solutions) < target_solutions):
    if enum_deadline is not None and time.time() >= enum_deadline:
        break
    prev_edge_sets = []
    for edge_set in enumerated_b:
        if edge_set and not any(edge_set == existing for existing in prev_edge_sets):
            prev_edge_sets.append(edge_set)
    next_sol = enumerate_next(prev_edge_sets, "nogood_hamming")
    if next_sol is None:
        break
    enumerated_b.append(next_sol["b_edges"])
    solutions.append(next_sol)
    finalize_solution(len(solutions) - 1, next_sol)
    if (not _unlimited) and len(solutions) >= target_solutions:
        break

if not solutions:
    raise RuntimeError("No feasible trees enumerated.")

def _best_selection_key(idx):
    sol = solutions[idx]
    fit_val = sol.get('residual_stage1') if score_bounded_migrations and args.migration_score_source == 'without-reward' else sol.get("fit")
    return (
        fit_val if fit_val is not None else float("inf"),
        -sol.get("root_degree", 0),
        idx
    )

best_idx = min(range(len(solutions)), key=_best_selection_key)

def _best_migration_key(idx):
    sol = solutions[idx]
    num_migrations = sol.get("num_migrations")
    num_comigrations = sol.get("num_comigrations")
    return (
        float("inf") if num_migrations is None else num_migrations,
        float("inf") if num_comigrations is None else num_comigrations,
        sol.get("fit", float("inf")),
        -sol.get("root_degree", 0),
        idx
    )

if migration_enabled:
    candidates = [i for i, sol in enumerate(solutions) if sol.get("migration_status") == "ok"]
    if not candidates:
        candidates = list(range(len(solutions)))
    best_mig_idx = min(candidates, key=_best_migration_key)
else:
    best_mig_idx = best_idx


# write all solutions and mark best-by-fit
with open(args.pathOutputFilePrefix + ".treeScores.tsv","w") as fs:
    fs.write("treeIndex\tfit\tfitWithoutReward\trootDegree\tnumMigrations\tnumComigrations\trefinedSolutionCount\trefinedEnumerationStatus\trefinedEnumerationExhaustive\tdagColoringCount\tisBestByFit\tisBestByMigration\tdagBeforePostprocess\tdagPostprocessStatus\tlabelingSource\trelabelChangeCount\tDAG\n")
    for i,sol in enumerate(solutions):
        B = sol["core"]["B"]
        rootDeg = sol.get("root_degree")
        if rootDeg is None:
            rootDeg = sum(1 for u in range(treeSize) if roundVar(B[0,u])==1)
        best_flag = 1 if i==best_idx else 0
        best_mig_flag = 1 if i==best_mig_idx else 0
        fit_val = sol.get("fit")
        fit_str = "NA" if fit_val is None else format(fit_val, '.17g')
        fit_no_reward_val = sol.get("residual_stage1")
        if fit_no_reward_val is None:
            try:
                core_i = sol["core"]
                fit_no_reward_val = core_i["residual_expr"].getValue() if core_i["residual_expr"] is not None else fit_val
            except Exception:
                fit_no_reward_val = fit_val
        fit_no_reward_str = "NA" if fit_no_reward_val is None else format(fit_no_reward_val, '.17g')
        mig_val = sol.get("num_migrations")
        mig_str = "NA" if mig_val is None else str(mig_val)
        comig_val = sol.get("num_comigrations")
        comig_str = "NA" if comig_val is None else str(comig_val)
        refined_count = sol.get("refined_solution_count", 0)
        refined_status = sol.get("refined_enumeration_status", "NA")
        refined_exhaustive = sol.get("refined_enumeration_exhaustive")
        refined_exhaustive_str = "NA" if refined_exhaustive is None else ("YES" if refined_exhaustive else "NO")
        dag_count = len(sol.get("migration_dag_solutions") or [])
        dag_before = sol.get("migration_dag_before_postprocess", "NA")
        dag_repair_status = sol.get("migration_dag_postprocess_status", "NA")
        labeling_source = sol.get("migration_labeling_source", "NA")
        relabel_change_count = sol.get("migration_relabel_change_count", "NA")
        dag_flag = migration_graph_dag_flag(sol)
        fs.write(f"{i}\t{fit_str}\t{fit_no_reward_str}\t{rootDeg}\t{mig_str}\t{comig_str}\t{refined_count}\t{refined_status}\t{refined_exhaustive_str}\t{dag_count}\t{best_flag}\t{best_mig_flag}\t{dag_before}\t{dag_repair_status}\t{labeling_source}\t{relabel_change_count}\t{dag_flag}\n")
print(f"[info] Wrote {len(solutions)} clonal tree solution(s); "
      f"best by DCF fit: {best_idx}; best by migration count: {best_mig_idx}.")

# Duplicate best solution with a convenient suffix
write_solution(best_mig_idx, solutions[best_mig_idx], tag="_best_migration")
 
