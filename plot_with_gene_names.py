#!/usr/bin/env python3
"""
Plot CITUPm tree file(s) with:
- Nodes colored by site from *.site_labels_sol{X}.tsv (matched to *.tree_{X}.txt).
- Optional override from embedded DAG_COLORING_INDEX blocks via --migration-dag-coloring-index.
- Optional batch rendering of all embedded CITUPm DAG colorings via --all-migration-dag-colorings.
- Top label with centered tree score/title; legend pinned top-right (inside label).
- 'Tree Score (no reward)' hidden if identical to 'Tree Score'.
- 'Samples:' removed (legend shows them).
 - Optional gene-role note nodes via --gene-info (Gene Symbol, Genome Location, Role in Cancer).
   Use --include-gene-other to also show genes without a role (labeled OTHER).
- Optional sample→site legend via --sample-sites (columns: sample, site).
- Optional site→color override via --site-colors (columns: site, color; hex colors recommended).
- NEW: Simple migration plotter for *.refined_edges*.tsv files; auto-finds matching refined_nodes and
  tree_{idx} files to map node ids to cluster ids and annotate sites.
- NEW: --with-site-transitions overlays site-to-site transitions (from matching migration_edges files,
  or embedded migration sections in the tree file) on main tree plots.
- NEW: --with_node_info adds "Node: <node_id>" to each non-root node in the main tree plots.

Usage:
  python plot_with_gene_names.py /path/to/score_mutation.tree_0.txt \
      --sample-sites /path/to/sample_site_map.tsv --gene-info /path/to/census.tsv
  python plot_with_gene_names.py /path/to/folder --sample-sites /path/to/sample_site_map.tsv
  python plot_with_gene_names.py /path/to/cn_neutr.refined_edges_sol1_best_migration.tsv
  python plot_with_gene_names.py /path/to/tree_3.txt --detopt /path/to/detopt_folder

Optional:
  --include-gene-other  # show genes lacking 'Role in Cancer' as OTHER
"""

import os
import sys
import glob
import re
import argparse
from collections import defaultdict
import pandas as pd
import graphviz


GENE_COLORS = ["red", "blue", "green", "orange", "purple", "brown", "magenta", "cyan"]
GENE_LINE_WRAP = 32  # max visible characters per gene line before wrapping

ROLE_TSG = "TSG"
ROLE_ONC = "Oncogene"
ROLE_OTHER = "Other"
ROLE_ORDER = [ROLE_TSG, ROLE_ONC, ROLE_OTHER]
ROLE_COLORS = {
    ROLE_TSG: "blue",
    ROLE_ONC: "red",
    ROLE_OTHER: "gray",
}
_GENE_SUFFIX_RE = re.compile(r"^(?P<base>.+)_(?P<idx>\d+)$")

# ---------------- helpers ----------------

def add_freq_lists(a, b):
    a_list = list(a) if a else []
    b_list = list(b) if b else []
    max_len = max(len(a_list), len(b_list))
    if max_len == 0:
        return []
    return [(a_list[i] if i < len(a_list) else 0.0) + (b_list[i] if i < len(b_list) else 0.0)
            for i in range(max_len)]


def parse_tree_file(filename):
    with open(filename, "r") as f:
        content = f.read().strip()
    sections = content.split("\n\n")

    tree_score = raw_tree_score = None
    lines = sections[0].splitlines()
    header_values = {}
    for i in range(0, len(lines) - 1, 2):
        key = lines[i].strip().upper()
        val = lines[i + 1].strip()
        if key:
            header_values[key] = val
    if lines and lines[0].strip().upper() == "TREE_SCORE":
        if len(lines) >= 2:
            tree_score = lines[1].strip()
        if len(lines) >= 4 and lines[2].strip().upper() == "TREE_SCORE_NO_REWARD":
            raw_tree_score = lines[3].strip()

    parent_child = []
    if len(sections) > 1:
        for line in sections[1].splitlines()[1:]:
            parts = line.split()
            if len(parts) >= 2:
                parent_child.append((parts[0], parts[1]))

    cluster_to_node, mutation_to_node, node_freqs, mutation_info = {}, {}, {}, {}

    if len(sections) > 2:
        for line in sections[2].splitlines()[1:]:
            parts = line.split()
            if len(parts) >= 2:
                cluster_to_node[parts[0]] = parts[1]

    if len(sections) > 3:
        for line in sections[3].splitlines()[1:]:
            parts = line.split()
            if len(parts) >= 2:
                mutation_to_node[parts[0]] = parts[1]

    if len(sections) > 4:
        for line in sections[4].splitlines()[1:]:
            parts = line.split()
            if len(parts) >= 3:
                node = parts[0]
                freqs = parts[1]
                try:
                    values = [float(x.strip()) for x in freqs.split(",") if x.strip()]
                    if values:
                        node_freqs[node] = values
                except Exception:
                    pass

    if len(sections) > 5:
        for line in sections[5].splitlines()[1:]:
            parts = line.split()
            if len(parts) >= 4:
                mutation_info[parts[0]] = (parts[1], parts[2], parts[3])

    origin_sites = {}
    dag_colorings = defaultdict(dict)
    dag_coloring_stats = {}
    migration_edges = []
    m_presence = defaultdict(dict)
    m_site_order = []
    seen_m_sites = set()

    for sec in sections[1:]:
        lines = [ln.strip() for ln in sec.splitlines() if ln.strip()]
        if not lines:
            continue

        header = lines[0].split()
        upper = [h.upper() for h in header]

        # Section style:
        # ORIGIN_SITE
        # NODE  SITE
        # ...
        if len(upper) == 1 and upper[0] == "ORIGIN_SITE":
            if len(lines) < 2:
                continue
            origin_header = lines[1].split()
            origin_upper = [h.upper() for h in origin_header]
            if "NODE" in origin_upper and "SITE" in origin_upper:
                node_idx = origin_upper.index("NODE")
                site_idx = origin_upper.index("SITE")
                need = max(node_idx, site_idx)
                for row in lines[2:]:
                    parts = row.split()
                    if len(parts) <= need:
                        continue
                    origin_sites[str(parts[node_idx])] = str(parts[site_idx])
            continue

        # Section style:
        # DAG_COLORING  SAME_ORIGIN_SCORE  MIGRATION_EDGE_COUNT  DAG
        # ...rows...
        if {"DAG_COLORING", "SAME_ORIGIN_SCORE", "MIGRATION_EDGE_COUNT", "DAG"}.issubset(set(upper)):
            idx_col = upper.index("DAG_COLORING")
            score_col = upper.index("SAME_ORIGIN_SCORE")
            edge_col = upper.index("MIGRATION_EDGE_COUNT")
            dag_col = upper.index("DAG")
            need = max(idx_col, score_col, edge_col, dag_col)
            for row in lines[1:]:
                parts = row.split()
                if len(parts) <= need:
                    continue
                try:
                    color_idx = int(parts[idx_col])
                except Exception:
                    continue
                dag_coloring_stats[color_idx] = {
                    "same_origin_score": parts[score_col],
                    "migration_edge_count": parts[edge_col],
                    "dag": parts[dag_col],
                }
            continue

        # Section style:
        # DAG_COLORING_INDEX NODE SITE
        # ...rows...
        if {"DAG_COLORING_INDEX", "NODE", "SITE"}.issubset(set(upper)):
            idx_col = upper.index("DAG_COLORING_INDEX")
            node_idx = upper.index("NODE")
            site_idx = upper.index("SITE")
            need = max(idx_col, node_idx, site_idx)
            for row in lines[1:]:
                parts = row.split()
                if len(parts) <= need:
                    continue
                try:
                    color_idx = int(parts[idx_col])
                except Exception:
                    continue
                dag_colorings[color_idx][str(parts[node_idx])] = str(parts[site_idx])
            continue

        # Section style:
        # REFINED_FROM_NODE REFINED_FROM_SITE REFINED_TO_NODE REFINED_TO_SITE EDGE_TYPE
        # ...rows...
        if {"REFINED_FROM_NODE", "REFINED_FROM_SITE", "REFINED_TO_NODE", "REFINED_TO_SITE"}.issubset(set(upper)):
            from_node_idx = upper.index("REFINED_FROM_NODE")
            from_site_idx = upper.index("REFINED_FROM_SITE")
            to_node_idx = upper.index("REFINED_TO_NODE")
            to_site_idx = upper.index("REFINED_TO_SITE")
            type_idx = upper.index("EDGE_TYPE") if "EDGE_TYPE" in upper else None
            need = max(from_node_idx, from_site_idx, to_node_idx, to_site_idx, type_idx or 0)
            for row in lines[1:]:
                parts = row.split()
                if len(parts) <= need:
                    continue
                etype = str(parts[type_idx]).strip().lower() if type_idx is not None else ""
                migration_edges.append({
                    "from_node": str(parts[from_node_idx]),
                    "from_site": str(parts[from_site_idx]),
                    "to_node": str(parts[to_node_idx]),
                    "to_site": str(parts[to_site_idx]),
                    "type": etype,
                })
            continue

        # Section style:
        # MIGRATION_FROM_NODE MIGRATION_FROM_SITE MIGRATION_TO_NODE MIGRATION_TO_SITE
        # ...rows...
        if {"MIGRATION_FROM_NODE", "MIGRATION_FROM_SITE", "MIGRATION_TO_NODE", "MIGRATION_TO_SITE"}.issubset(set(upper)):
            from_node_idx = upper.index("MIGRATION_FROM_NODE")
            from_site_idx = upper.index("MIGRATION_FROM_SITE")
            to_node_idx = upper.index("MIGRATION_TO_NODE")
            to_site_idx = upper.index("MIGRATION_TO_SITE")
            need = max(from_node_idx, from_site_idx, to_node_idx, to_site_idx)
            for row in lines[1:]:
                parts = row.split()
                if len(parts) <= need:
                    continue
                migration_edges.append({
                    "from_node": str(parts[from_node_idx]),
                    "from_site": str(parts[from_site_idx]),
                    "to_node": str(parts[to_node_idx]),
                    "to_site": str(parts[to_site_idx]),
                    "type": "migration",
                })
            continue

        # Section style:
        # REFINED_NODE SITE PRESENCE_SIGNAL OBSERVED MU D IS_ORIGIN
        # ...rows...
        if {"REFINED_NODE", "SITE", "MU"}.issubset(set(upper)):
            node_idx = upper.index("REFINED_NODE")
            site_idx = upper.index("SITE")
            m_idx = upper.index("MU")
            need = max(node_idx, site_idx, m_idx)
            for row in lines[1:]:
                parts = row.split()
                if len(parts) <= need:
                    continue
                node = str(parts[node_idx])
                site = str(parts[site_idx])
                try:
                    m_val = float(parts[m_idx])
                    m_bit = 1 if m_val > 0 else 0
                except Exception:
                    m_bit = 0
                m_presence[node][site] = m_bit
                if site and site not in seen_m_sites:
                    seen_m_sites.add(site)
                    m_site_order.append(site)
            continue

        # Section style:
        # MU_ASSIGNMENT SITE MU
        # ...rows...
        if {"MU_ASSIGNMENT", "SITE", "MU"}.issubset(set(upper)):
            node_idx = upper.index("MU_ASSIGNMENT")
            site_idx = upper.index("SITE")
            m_idx = upper.index("MU")
            need = max(node_idx, site_idx, m_idx)
            for row in lines[1:]:
                parts = row.split()
                if len(parts) <= need:
                    continue
                node = str(parts[node_idx])
                site = str(parts[site_idx])
                try:
                    m_val = float(parts[m_idx])
                    m_bit = 1 if m_val > 0 else 0
                except Exception:
                    m_bit = 0
                m_presence[node][site] = m_bit
                if site and site not in seen_m_sites:
                    seen_m_sites.add(site)
                    m_site_order.append(site)

    return {
        "tree_score": tree_score,
        "raw_tree_score": raw_tree_score,
        "header_values": header_values,
        "parent_child": parent_child,
        "cluster_to_node": cluster_to_node,
        "mutation_to_node": mutation_to_node,
        "node_freqs": node_freqs,
        "mutation_info": mutation_info,
        "origin_sites": origin_sites,
        "dag_colorings": {idx: dict(site_map) for idx, site_map in dag_colorings.items()},
        "dag_coloring_stats": dict(dag_coloring_stats),
        "migration_edges": migration_edges,
        "m_presence": {node: dict(site_map) for node, site_map in m_presence.items()},
        "m_site_order": list(m_site_order),
    }


def choose_origin_sites(data, dag_coloring_index=None):
    if dag_coloring_index is None:
        return dict(data.get("origin_sites") or {}), None
    dag_colorings = data.get("dag_colorings") or {}
    chosen = dag_colorings.get(int(dag_coloring_index))
    if not chosen:
        return None, None
    return dict(chosen), int(dag_coloring_index)


def derive_migration_transitions(origin_sites, m_presence, m_site_order=None):
    transitions = defaultdict(list)
    transition_sites = set()
    ordered_sites = list(m_site_order or [])
    ordered_site_set = set(ordered_sites)

    for site_map in (m_presence or {}).values():
        for site in site_map:
            if site and site not in ordered_site_set:
                ordered_site_set.add(site)
                ordered_sites.append(site)

    for node_id, site_map in (m_presence or {}).items():
        node_key = str(node_id)
        if node_key == "0":
            continue
        origin = (origin_sites or {}).get(node_key)
        if not origin:
            continue
        for site in ordered_sites:
            if site == origin or not site_map.get(site, 0):
                continue
            transitions[node_key].append({
                "from_site": origin,
                "to_site": site,
                "type": "migration",
            })
            transition_sites.add(origin)
            transition_sites.add(site)
    return transitions, transition_sites


def effective_presence_bits(node_id, site_map, ordered_sites, origin_sites=None):
    node_key = str(node_id)
    if node_key == "0" and origin_sites and origin_sites.get(node_key):
        root_site = origin_sites.get(node_key)
        return [1 if site == root_site else 0 for site in ordered_sites]
    return [1 if site_map.get(site, 0) else 0 for site in ordered_sites]


def compute_node_stats(data):
    # children
    children = defaultdict(list)
    for p, c in data["parent_child"]:
        children[str(p)].append(str(c))

    nodes = set()
    for p, c in data["parent_child"]:
        nodes.add(str(p)); nodes.add(str(c))
    nodes.update(str(n) for n in data["node_freqs"].keys())
    nodes.update(str(n) for n in data["cluster_to_node"].values())
    nodes.update(str(n) for n in data["mutation_to_node"].values())

    # average observed per node
    obs_sum, obs_cnt = {}, {}
    for mut, (obs, _, _) in data["mutation_info"].items():
        n = data["mutation_to_node"].get(mut)
        if not n:
            continue
        n = str(n)
        try:
            vals = [float(x.strip()) for x in obs.split(",") if x.strip()]
        except Exception:
            continue
        if n not in obs_sum:
            obs_sum[n] = [0.0] * len(vals)
            obs_cnt[n] = 0
        elif len(obs_sum[n]) != len(vals):
            continue
        for i, v in enumerate(vals):
            obs_sum[n][i] += v
        obs_cnt[n] += 1

    avg_obs = {n: (tuple(v / obs_cnt[n] for v in obs_sum[n]) if n in obs_cnt and obs_cnt[n] else None) for n in nodes}
    node_to_cluster = {node: cl for cl, node in data["cluster_to_node"].items()}

    # lineage totals
    total_freqs = {}
    def acc(n):
        if n in total_freqs:
            return total_freqs[n]
        base = data["node_freqs"].get(n)
        total = list(base) if base else []
        for ch in children.get(n, []):
            total = add_freq_lists(total, acc(ch))
        total_freqs[n] = total
        return total
    for n in sorted_nodes_by_cluster(nodes, node_to_cluster):
        acc(n)

    stats = {}
    for n in sorted_nodes_by_cluster(nodes, node_to_cluster):
        stats[n] = {
            "node": n,
            "base": data["node_freqs"].get(n),
            "lineage": total_freqs.get(n),
            "observed": avg_obs.get(n),
        }
    return stats


def parse_gene_info(file_path, include_no_role=False):
    try:
        df = pd.read_csv(file_path, sep="\t")
    except Exception as e:
        sys.stderr.write(f"Error reading gene info file {file_path}: {e}\n")
        return {}
    for col in ["Gene Symbol", "Genome Location", "Role in Cancer"]:
        if col not in df.columns:
            sys.stderr.write(f"Gene info file must contain '{col}' column.\n")
            return {}
    synonym_col = "Synonyms" if "Synonyms" in df.columns else None
    entries = []
    for _, row in df.iterrows():
        gene = str(row["Gene Symbol"]).strip()
        role = str(row["Role in Cancer"]).strip()
        role_missing = role.lower() in ["", "nan"]
        if role_missing and not include_no_role:
            continue
        loc = str(row["Genome Location"]).strip()
        try:
            chrom, coords = loc.split(":", 1)
            start_str, end_str = coords.split("-")
            entry = {
                "base": gene,
                "role": role if not role_missing else "",
                "chrom": chrom,
                "start": int(start_str),
                "end": int(end_str),
            }
            if synonym_col:
                syn_field = "" if pd.isna(row.get(synonym_col)) else str(row[synonym_col])
                synonyms = [s.strip() for s in syn_field.split(",") if s.strip()]
                entry["synonyms"] = synonyms
            entries.append(entry)
        except Exception:
            sys.stderr.write(f"Bad Genome Location '{loc}' for {gene}\n")
    if not entries:
        return {}
    counts = defaultdict(int)
    for entry in entries:
        counts[entry["base"]] += 1
    out = {}
    seen = defaultdict(int)
    for entry in entries:
        base = entry["base"]
        if counts[base] > 1:
            seen[base] += 1
            label = f"{base}_{seen[base]}"
        else:
            label = base
        out[label] = entry
    return out


def classify_gene_roles(role_text, include_other=False):
    roles = set()
    if role_text is None:
        role_text = ""
    role_str = str(role_text).strip()
    role_lower = role_str.lower()
    if not role_lower or role_lower == "nan":
        if include_other:
            roles.add(ROLE_OTHER)
        return roles
    has_fusion = "fusion" in role_lower
    if "tsg" in role_lower:
        roles.add(ROLE_TSG)
    if "oncogene" in role_lower:
        roles.add(ROLE_ONC)
    if not roles and include_other and not has_fusion:
        roles.add(ROLE_OTHER)
    return roles


def build_gene_role_map(gene_info):
    role_map = defaultdict(set)
    if not gene_info:
        return {}
    for label, info in gene_info.items():
        base = info.get("base", label)
        roles = classify_gene_roles(info.get("role", ""), include_other=False)
        if not roles:
            continue
        role_map[base.lower()].update(roles)
        for syn in info.get("synonyms", []) if isinstance(info, dict) else []:
            role_map[syn.lower()].update(roles)
    return role_map


def gene_base_name(label):
    if not label:
        return label
    m = _GENE_SUFFIX_RE.match(str(label))
    return m.group("base") if m else str(label)


def group_genes_by_role(genes, gene_role_map, include_other=False):
    grouped = defaultdict(set)
    role_map = gene_role_map or {}
    for g in genes:
        g_lower = g.lower()
        roles = role_map.get(g_lower, set())
        if not roles:
            base = gene_base_name(g).lower()
            roles = role_map.get(base, set())
        if not roles and include_other:
            roles = {ROLE_OTHER}
        for role in roles:
            grouped[role].add(g)
    return grouped


def parse_mutation_coordinate(mutation_id):
    """
    Mutation ids sometimes include the ref/alt allele (e.g., 'chr1:12345:A').
    We only need the chromosome and position for gene matching, so accept ids
    with 2 or 3 colon-separated parts.
    """
    try:
        parts = mutation_id.split(":")
        if len(parts) < 2:
            return None, None
        chrom, pos = parts[0], parts[1]
        if chrom.lower().startswith("chr"):
            chrom = chrom[3:]
        return chrom, int(pos)
    except Exception:
        return None, None


def build_gene_note_html(role_to_genes):
    rows = []
    for i, role in enumerate(sorted(role_to_genes)):
        genes = ", ".join(sorted(role_to_genes[role]))
        rows.append(f'<TR><TD><FONT COLOR="{GENE_COLORS[i % len(GENE_COLORS)]}"><B>{role.upper()}</B></FONT>: {genes}</TD></TR>')
    return "<<TABLE BORDER='0' CELLBORDER='0' CELLSPACING='4'>" + "".join(rows) + "</TABLE>>"


def wrap_role_gene_lines(role, genes, color, max_chars=GENE_LINE_WRAP):
    """
    Wrap gene names for a role into multiple lines (no role label included) to avoid very wide nodes.
    Returns list of HTML-safe strings (without surrounding <TR>/<TD>).
    """
    def visible_len(s):
        return len(s.replace("&nbsp;", " "))

    lines = []
    cur = ""
    for g in sorted(genes):
        token = g if not cur else ", " + g
        if visible_len(cur + token) > max_chars and cur:
            lines.append(cur)
            cur = g
        else:
            cur += token
    lines.append(cur)
    return lines


def allowed_gene_set(gene_info):
    if not gene_info:
        return None
    allowed = set()
    for gene, info in gene_info.items():
        base = info.get("base", gene)
        allowed.add(gene.lower())
        allowed.add(base.lower())
        for syn in info.get("synonyms", []) if isinstance(info, dict) else []:
            allowed.add(syn.lower())
    return allowed


def format_states(states):
    if not states:
        return ""
    cleaned = []
    for s in states:
        if s is None:
            continue
        ss = str(s).strip()
        if not ss:
            continue
        if ss.startswith("(") and ss.endswith(")"):
            cleaned.append(ss)
        else:
            cleaned.append(f"({ss})")
    return "; ".join(cleaned)


# --------------- sites & legend ---------------

SITE_COLORS = [
    "#add8e6", "#98fb98", "#f0e68c", "#e5a79e", "#a5b1eb",
    "#ccf0ad", "#f4b6e5", "#92f2ed", "#e8d2a8", "#c7b0ed",
    "#b8f2bc", "#ef94a8", "#9dc9f3", "#dfebb2", "#ebbaf0",
    "#98eccc", "#f1ba9f", "#aba9f4", "#c9edbd", "#e99bc6",
    "#a2e2ee", "#f3e9aa", "#dab4f6", "#9ee5b2", "#eba5a6",
    "#adc1f0", "#dbf4b6", "#f292e6", "#a8e8dd", "#edd1b0",
    "#c6b8f2", "#99ef94", "#f39dbb", "#b2d6eb", "#ecf0ba",
    "#da98ec", "#9ff1c7", "#f4b8a9", "#bdc1ed", "#b9e99b",
    "#eea2d6", "#aaf0f3", "#f6e5b4", "#be9ee5", "#a5ebb0",
    "#f0adb6", "#b6d0f4", "#d7f292", "#e7a8e8", "#b0eddb",
    "#f2cfb8", "#9f94ef", "#adf39d", "#ebb2cd", "#bae3f0",
    "#ece898", "#d59ff1", "#a9f4c4", "#edc0bd", "#9bace9",
    "#c9eea2", "#f3aae5", "#b4f6ef", "#e5ca9e", "#bca5eb",
    "#adf0af", "#f4b6c6", "#92c7f2", "#dee8a8", "#e5b0ed",
    "#b8f2d9", "#efae94", "#9d9ff3", "#c4ebb2", "#f0bada",
    "#98e3ec", "#f1e29f", "#d1a9f4", "#bdedc8", "#e99b9f",
    "#a2bdee", "#d9f3aa", "#f6b4f1", "#9ee5d5", "#ebc7a5",
    "#baadf0", "#bcf4b6", "#f292b7", "#a8d4e8", "#ecedb0",
    "#e2b8f2", "#94efbd", "#f3aa9d", "#b2baeb", "#d1f0ba",
    "#ec98d5", "#9ff1f0", "#f4dda9", "#d0bded", "#9be9a3",
]

def parse_site_labels(path):
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as e:
        sys.stderr.write(f"Error reading site labels file {path}: {e}\n")
        return {}
    if not {"node", "site"}.issubset(df.columns):
        sys.stderr.write(f"Site labels file {path} must have columns 'node' and 'site'.\n")
        return {}
    out = {}
    for _, row in df.iterrows():
        out[str(row["node"]).strip()] = "" if pd.isna(row["site"]) else str(row["site"]).strip()
    return out


def parse_sample_sites(path):
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as e:
        sys.stderr.write(f"Error reading sample-sites file {path}: {e}\n")
        return {}
    if not {"sample", "site"}.issubset(df.columns):
        sys.stderr.write(f"Sample-sites file {path} must have columns 'sample' and 'site'.\n")
        return {}
    out = {}
    for _, row in df.iterrows():
        s = str(row["sample"]).strip()
        if s:
            out[s] = "" if pd.isna(row["site"]) else str(row["site"]).strip()
    return out


def parse_site_colors(path):
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as e:
        sys.stderr.write(f"Error reading site-colors file {path}: {e}\n")
        return {}
    if not {"site", "color"}.issubset(df.columns):
        sys.stderr.write(f"Site-colors file {path} must have columns 'site' and 'color'.\n")
        return {}
    out = {}
    hex_re = re.compile(r"^#[0-9A-Fa-f]{6}$")
    name_re = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
    for _, row in df.iterrows():
        site = "" if pd.isna(row["site"]) else str(row["site"]).strip()
        color = "" if pd.isna(row["color"]) else str(row["color"]).strip()
        if not site or not color:
            continue
        if not (hex_re.match(color) or name_re.match(color)):
            sys.stderr.write(
                f"Warning: skipping color '{color}' for site '{site}'. "
                "Use #RRGGBB hex or a Graphviz color name.\n"
            )
            continue
        out[site] = color
    return out


def natural_sort_key(value):
    parts = re.split(r"(\d+)", str(value))
    return [int(p) if p.isdigit() else p.lower() for p in parts]


def ordered_sample_site_items(sample_site_map):
    if not sample_site_map:
        return []
    return sorted(sample_site_map.items(), key=lambda x: natural_sort_key(x[0]))


def sample_site_items_for_width(sample_site_map, width):
    if width <= 0:
        return []
    if sample_site_map and all(f"sample{i}" in sample_site_map for i in range(width)):
        return [(f"sample{i}", sample_site_map.get(f"sample{i}", "")) for i in range(width)]

    items = ordered_sample_site_items(sample_site_map)
    if len(items) >= width:
        return items[:width]

    out = list(items)
    for i in range(len(out), width):
        out.append((f"sample{i}", ""))
    return out


def find_sample_sites_file_for_tree(tree_path):
    patient = os.path.basename(tree_path).split(".", 1)[0]
    if not patient:
        return None

    cur = os.path.abspath(os.path.dirname(tree_path))
    ancestors = []
    while True:
        ancestors.append(cur)
        parent = os.path.dirname(cur)
        if parent == cur:
            break
        cur = parent

    seen = set()
    for root in ancestors:
        candidate_dirs = [
            root,
            os.path.join(root, "mig8_sample_sites_machina_style"),
            os.path.join(root, "tracerX_data", "zenodo_7649257", "mig8_sample_sites_machina_style"),
        ]
        for d in candidate_dirs:
            path = os.path.join(d, f"{patient}.sites.tsv")
            if path in seen:
                continue
            seen.add(path)
            if os.path.exists(path):
                return path
    return None


def make_site_color_map(*site_sets, site_order=None, override_colors=None):
    sites = set()
    for site in site_order or []:
        if site:
            sites.add(site)
    for sset in site_sets:
        for site in sset:
            if site:
                sites.add(site)
    uniques = sorted(sites, key=natural_sort_key)
    color_map = {site: SITE_COLORS[i % len(SITE_COLORS)] for i, site in enumerate(uniques)}
    for site, color in (override_colors or {}).items():
        if site in color_map:
            color_map[site] = color
    return color_map


def sort_value_key(value):
    parts = re.split(r"(\d+)", str(value))
    out = []
    for part in parts:
        if part == "":
            continue
        out.append((0, int(part)) if part.isdigit() else (1, part.lower()))
    return tuple(out)


def node_cluster_sort_key(node_id, node_to_cluster=None):
    node = str(node_id)
    if node == "0":
        return (-1, (), sort_value_key(node))
    cluster = (node_to_cluster or {}).get(node)
    if cluster is not None:
        return (0, sort_value_key(cluster), sort_value_key(node))
    return (1, sort_value_key(node), sort_value_key(node))


def sorted_nodes_by_cluster(nodes, node_to_cluster=None):
    return sorted((str(n) for n in nodes), key=lambda n: node_cluster_sort_key(n, node_to_cluster))


def sorted_parent_child_edges(parent_child, node_to_cluster=None):
    children = defaultdict(list)
    parents, child_nodes = set(), set()
    for p, c in parent_child:
        p, c = str(p), str(c)
        children[p].append(c)
        parents.add(p)
        child_nodes.add(c)

    roots = sorted(parents - child_nodes, key=lambda n: node_cluster_sort_key(n, node_to_cluster))
    if "0" in parents and "0" not in roots:
        roots.insert(0, "0")

    emitted = []
    seen_edges = set()

    def emit_from(parent):
        for child in sorted(children.get(parent, []), key=lambda n: node_cluster_sort_key(n, node_to_cluster)):
            edge = (parent, child)
            if edge in seen_edges:
                continue
            seen_edges.add(edge)
            emitted.append(edge)
            emit_from(child)

    for root in roots:
        emit_from(root)

    for p, c in parent_child:
        edge = (str(p), str(c))
        if edge not in seen_edges:
            seen_edges.add(edge)
            emitted.append(edge)

    return emitted


def parse_tree_filename(tree_path):
    m = re.match(r"(?P<prefix>.+)\.tree_(?P<idx>\d+)(?P<suffix>.*)\.txt$", os.path.basename(tree_path))
    if not m:
        return None, None, None
    return m.group("prefix"), m.group("idx"), m.group("suffix")


def extract_solution_index_from_tree(tree_path):
    _, idx, _ = parse_tree_filename(tree_path)
    return idx


def is_best_migration_tree(tree_path):
    _, _, suffix = parse_tree_filename(tree_path)
    return bool(suffix and "best_migration" in suffix)


def discover_tree_files(path, include_best_migration=False):
    tree_files = sorted(glob.glob(os.path.join(path, "*.tree_*.txt")))
    if not include_best_migration:
        tree_files = [p for p in tree_files if not is_best_migration_tree(p)]
    return tree_files


def find_migration_edges_for_tree(tree_path):
    """
    Try to locate the migration_edges file corresponding to a given tree_{idx}*.txt.
    """
    prefix, idx, suffix = parse_tree_filename(tree_path)
    if not prefix:
        return None
    d = os.path.dirname(tree_path)
    candidates = [
        os.path.join(d, f"{prefix}.refined_edges_sol{idx}{suffix}.tsv"),
    ]
    # Backward compatibility for older mig8 outputs: tree_0.txt used unsuffixed
    # CITUPm sidecars, while tree_0_best_migration.txt already used _sol0.
    if idx == "0" and suffix == "":
        candidates.append(os.path.join(d, f"{prefix}.refined_edges.tsv"))
    for direct in candidates:
        if os.path.exists(direct):
            return direct
    hits = sorted(glob.glob(os.path.join(d, f"{prefix}.refined_edges_sol{idx}{suffix}*.tsv")))
    return hits[0] if hits else None


def find_site_labels_file_for_tree(tree_path):
    prefix, idx, suffix = parse_tree_filename(tree_path)
    if idx is None:
        return None
    d = os.path.dirname(tree_path)
    exact_hits = []
    if prefix:
        exact_hits = glob.glob(os.path.join(d, f"{prefix}.site_labels_sol{idx}{suffix}.tsv"))
    hits = exact_hits or glob.glob(os.path.join(d, f"*.site_labels_sol{idx}{suffix}*.tsv"))
    if not hits and idx == "0" and suffix == "":
        hits = glob.glob(os.path.join(d, "*.site_labels.tsv"))
    return hits[0] if hits else None


# ---------------- refined migration helpers ----------------

def _safe_float(val):
    try:
        if pd.isna(val):
            return None
        return float(val)
    except Exception:
        return None


def _clean_node_id(val):
    if pd.isna(val):
        return None
    txt = str(val).strip()
    try:
        num = float(txt)
        if num.is_integer():
            return str(int(num))
    except Exception:
        pass
    return txt


def parse_migration_edges(path):
    required = {"from_node", "from_site", "to_node", "to_site", "type"}
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as e:
        sys.stderr.write(f"Error reading refined edges file {path}: {e}\n")
        return []
    if not required.issubset(df.columns):
        sys.stderr.write(f"CITUPm edges file {path} must have columns: {', '.join(sorted(required))}\n")
        return []
    edges = []
    for _, row in df.iterrows():
        from_node = _clean_node_id(row["from_node"])
        to_node = _clean_node_id(row["to_node"])
        if from_node is None or to_node is None:
            continue
        edges.append({
            "from_node": from_node,
            "from_site": "" if pd.isna(row["from_site"]) else str(row["from_site"]).strip(),
            "to_node": to_node,
            "to_site": "" if pd.isna(row["to_site"]) else str(row["to_site"]).strip(),
            "type": str(row["type"]).strip().lower() if not pd.isna(row["type"]) else "",
        })
    return edges


def parse_migration_nodes(path):
    required = {"node", "site"}
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as e:
        sys.stderr.write(f"Error reading refined nodes file {path}: {e}\n")
        return {}
    if not required.issubset(df.columns):
        sys.stderr.write(f"CITUPm nodes file {path} must have columns at least: {', '.join(sorted(required))}\n")
        return {}
    out = defaultdict(list)
    for _, row in df.iterrows():
        node = _clean_node_id(row["node"])
        if node is None:
            continue
        site = "" if pd.isna(row["site"]) else str(row["site"]).strip()
        origin_val = _safe_float(row["is_origin"]) if "is_origin" in df.columns else None
        out[node].append({
            "site": site,
            "presence_signal": _safe_float(row["presence_signal"]) if "presence_signal" in df.columns else None,
            "observed": _safe_float(row["observed"]) if "observed" in df.columns else None,
            "is_origin": bool(origin_val) if origin_val is not None else False,
        })
    return out


def site_order_from_migration_nodes(path):
    """
    Return the first-seen order of non-empty site names in a migration_nodes file.
    """
    try:
        df = pd.read_csv(path, sep="\t", usecols=["site"])
    except Exception:
        return []
    seen, order = set(), []
    for v in df["site"]:
        if pd.isna(v):
            continue
        s = str(v).strip()
        if s and s not in seen:
            seen.add(s)
            order.append(s)
    return order


# ---------------- detopt helpers ----------------

def _format_number(val):
    try:
        num = float(val)
    except Exception:
        return str(val)
    if num.is_integer():
        return str(int(num))
    return f"{num:g}"


def _parse_copy_pair(val):
    if val is None or pd.isna(val):
        return None
    txt = str(val).strip()
    if not txt:
        return None
    parts = [p.strip() for p in txt.split(",")]
    if len(parts) < 2:
        return None
    try:
        return float(parts[0]), float(parts[1])
    except Exception:
        return None


def load_detopt_copy_numbers(path):
    """
    Returns dict: mut_id -> {cluster: (allele_a, allele_b)}
    """
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as e:
        sys.stderr.write(f"Error reading detopt copy numbers file {path}: {e}\n")
        return {}
    if "mut_id" not in df.columns:
        sys.stderr.write(f"Detopt copy numbers file {path} missing required 'mut_id' column.\n")
        return {}
    out = {}
    cols = [c for c in df.columns if c != "mut_id"]
    for _, row in df.iterrows():
        mut_id = "" if pd.isna(row.get("mut_id")) else str(row["mut_id"]).strip()
        if not mut_id:
            continue
        row_map = {}
        for col in cols:
            pair = _parse_copy_pair(row.get(col))
            if pair is None:
                continue
            row_map[str(col).strip()] = pair
        if row_map:
            out[mut_id] = row_map
    return out


def _state_numbers_from_text(state):
    if state is None:
        return None
    nums = re.findall(r"-?\d+(?:\.\d+)?", str(state))
    return nums if nums else None


def _format_state_from_numbers(nums, allele_marks=None):
    if not nums:
        return None
    out = [str(n) for n in nums]
    if allele_marks:
        if len(out) >= 1 and allele_marks[0]:
            out[0] = f"{out[0]}*"
        if len(out) >= 2 and allele_marks[1]:
            out[1] = f"{out[1]}*"
    return "(" + ", ".join(out) + ")"


def _allele_marks_for_mutation(mut_copy_map, mut_key, cluster):
    if not mut_copy_map or not mut_key:
        return None
    pair = mut_copy_map.get(mut_key, {}).get(str(cluster))
    if not pair:
        return None
    return (pair[0] > 0, pair[1] > 0)


def _format_state_with_allele(state, mut_key, cluster, mut_copy_map, tot_copy_map, apply_allele_marks=True):
    allele_marks = _allele_marks_for_mutation(mut_copy_map, mut_key, cluster) if apply_allele_marks else None
    nums = _state_numbers_from_text(state) if state else None
    if not nums and tot_copy_map and mut_key:
        totals = tot_copy_map.get(mut_key, {}).get(str(cluster))
        if totals:
            nums = [_format_number(v) for v in totals]
    if nums:
        return _format_state_from_numbers(nums, allele_marks)
    return state

def _mutation_key_from_parts(chr_val, start, end, ref, alt):
    return f"chr{chr_val}:{start}-{end}_{ref}>{alt}"


def load_detopt_multianno(path):
    """
    Returns dict: mut_key -> set(genes)
    """
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as e:
        sys.stderr.write(f"Error reading detopt multianno file {path}: {e}\n")
        return {}
    required = {"Chr", "Start", "End", "Ref", "Alt", "Gene.refGene"}
    if not required.issubset(df.columns):
        sys.stderr.write(f"Detopt multianno file {path} missing required columns.\n")
        return {}
    gene_locations = defaultdict(set)
    for _, row in df.iterrows():
        try:
            chrom = row["Chr"]
            start = row["Start"]
            end = row["End"]
        except Exception:
            continue
        gene_field = "" if pd.isna(row["Gene.refGene"]) else str(row["Gene.refGene"]).strip()
        for chunk in re.split(r"[;,]", gene_field):
            g = chunk.strip()
            if g and g != ".":
                gene_locations[g].add((chrom, start, end))

    gene_loc_labels = {}
    for gene, locs in gene_locations.items():
        if len(locs) <= 1:
            for loc in locs:
                gene_loc_labels[(gene, loc)] = gene
            continue
        def sort_key(loc):
            chrom, start, end = loc
            chrom_txt = str(chrom)
            if chrom_txt.lower().startswith("chr"):
                chrom_txt = chrom_txt[3:]
            chrom_key = (0, int(chrom_txt)) if chrom_txt.isdigit() else (1, chrom_txt)
            try:
                start_key = int(start)
            except Exception:
                start_key = str(start)
            try:
                end_key = int(end)
            except Exception:
                end_key = str(end)
            return (chrom_key, start_key, end_key)
        for idx, loc in enumerate(sorted(locs, key=sort_key), start=1):
            gene_loc_labels[(gene, loc)] = f"{gene}_{idx}"

    mut_to_genes = {}
    for _, row in df.iterrows():
        try:
            key = _mutation_key_from_parts(row["Chr"], row["Start"], row["End"], row["Ref"], row["Alt"])
            chrom = row["Chr"]
            start = row["Start"]
            end = row["End"]
        except Exception:
            continue
        gene_field = "" if pd.isna(row["Gene.refGene"]) else str(row["Gene.refGene"]).strip()
        genes = set()
        loc = (chrom, start, end)
        for chunk in re.split(r"[;,]", gene_field):
            g = chunk.strip()
            if g and g != ".":
                genes.add(gene_loc_labels.get((g, loc), g))
        if not genes:
            continue
        mut_to_genes.setdefault(key, set()).update(genes)
    return mut_to_genes


def load_detopt_assignments(path):
    """
    Returns list of dicts per mutation row.
    """
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as e:
        sys.stderr.write(f"Error reading detopt assignments file {path}: {e}\n")
        return []
    cols = set(df.columns)
    has_snv = "snv_assignment" in cols
    cna_assign_cols = [c for c in cols if c.startswith("cna_assignment")]
    assignments = []
    for _, row in df.iterrows():
        entry = {
            "mut_index": "" if pd.isna(row.get("mut_index")) else str(row["mut_index"]).strip(),
            "snv_assignment": row["snv_assignment"] if has_snv and not pd.isna(row["snv_assignment"]) else None,
            "cna": [],
        }
        for cna_col in cna_assign_cols:
            suffix = cna_col.replace("cna_assignment", "")
            cluster_val = row[cna_col]
            if pd.isna(cluster_val):
                continue
            suffix_clean = suffix.strip("_")
            # state column can be "cna_state0" or "cna_state_0" etc.
            candidate_state_cols = [
                f"cna_state{suffix}",
                f"cna_state_{suffix_clean}",
            ]
            state_val = None
            for sc in candidate_state_cols:
                if sc in df.columns:
                    state_val = row[sc]
                    break
            entry["cna"].append((suffix_clean, cluster_val, state_val))
        assignments.append(entry)
    return assignments


def build_detopt_cluster_genes(assignments_path, multianno_path, allowed_genes=None,
                               mut_copy_path=None, tot_copy_path=None):
    mut_to_genes = load_detopt_multianno(multianno_path)
    rows = load_detopt_assignments(assignments_path)
    mut_copy_map = load_detopt_copy_numbers(mut_copy_path) if mut_copy_path else {}
    tot_copy_map = load_detopt_copy_numbers(tot_copy_path) if tot_copy_path else {}
    cluster_map = defaultdict(lambda: defaultdict(lambda: {"snv": False, "cna_states": {}}))
    cluster_mut_counts = defaultdict(int)
    for row in rows:
        mut_key = row.get("mut_index")
        genes = mut_to_genes.get(mut_key, None)
        if not genes:
            # fallback to mut id as pseudo gene
            genes = {mut_key}
        if allowed_genes is None:
            genes_filtered = genes
        else:
            genes_filtered = {
                g for g in genes
                if g and (g.lower() in allowed_genes or gene_base_name(g).lower() in allowed_genes)
            }
        # snv assignment
        snv_cl = row.get("snv_assignment")
        snv_cluster = None
        if snv_cl is not None:
            try:
                snv_cl = str(int(snv_cl))
            except Exception:
                snv_cl = str(snv_cl).strip()
            snv_cluster = snv_cl
            cluster_mut_counts[snv_cl] += 1
            if genes_filtered:
                for g in genes_filtered:
                    info = cluster_map[snv_cl][g]
                    info["snv"] = True
                    snv_state = _format_state_with_allele(
                        None,
                        mut_key,
                        snv_cl,
                        mut_copy_map,
                        tot_copy_map,
                        apply_allele_marks=True,
                    )
                    if snv_state and "snv" not in info["cna_states"]:
                        info["cna_states"]["snv"] = snv_state
        # cna assignments
        for suffix, cluster_val, state_val in row.get("cna", []):
            try:
                c_cl = str(int(cluster_val))
            except Exception:
                c_cl = str(cluster_val).strip()
            cluster_mut_counts[c_cl] += 1
            state_raw = None if pd.isna(state_val) else str(state_val).strip()
            apply_marks = snv_cluster is not None and c_cl == snv_cluster
            state = _format_state_with_allele(
                state_raw,
                mut_key,
                c_cl,
                mut_copy_map,
                tot_copy_map,
                apply_allele_marks=apply_marks,
            )
            if genes_filtered:
                for g in genes_filtered:
                    info = cluster_map[c_cl][g]
                    if state:
                        info["cna_states"][suffix] = state
                        if "snv" in info["cna_states"]:
                            del info["cna_states"]["snv"]
    output = {}
    for cl, genes in cluster_map.items():
        output[cl] = {g: {"snv": info["snv"], "cna_states": dict(info["cna_states"])} for g, info in genes.items()}
    return output, cluster_mut_counts


def find_detopt_files(detopt_dir, tree_path):
    """
    Return (assignments_path, multianno_path, mut_copy_path, tot_copy_path) matching the tree index, if found.
    """
    if not detopt_dir:
        return None, None, None, None
    base = os.path.basename(tree_path)
    m = re.match(r"(?P<prefix>.+)\.tree_(?P<idx>\d+)(?P<suffix>.*)\.txt$", base)
    prefix = m.group("prefix") if m else None
    idx = m.group("idx") if m else extract_solution_index_from_tree(tree_path)
    suffix = m.group("suffix") if m else ""
    if idx is None:
        return None, None, None, None

    def collect_hits(patterns):
        hits = []
        for pat in patterns:
            hits.extend(glob.glob(os.path.join(detopt_dir, pat)))
        return sorted(set(hits))

    assign_hits = []
    anno_hits = []
    mut_copy_hits = []
    tot_copy_hits = []
    if prefix:
        assign_hits = collect_hits([
            f"{prefix}.tree_{idx}{suffix}_assignments.tsv",
            f"{prefix}.tree_{idx}{suffix}*_assignments.tsv",
            f"{prefix}.tree_{idx}*_assignments.tsv",
        ])
        anno_hits = collect_hits([
            f"{prefix}.tree_{idx}{suffix}.hg38_multianno.txt",
            f"{prefix}.tree_{idx}{suffix}*.hg38_multianno.txt",
            f"{prefix}.tree_{idx}*.hg38_multianno.txt",
        ])
        mut_copy_hits = collect_hits([
            f"{prefix}.tree_{idx}{suffix}_mut_copy_numbers.tsv",
            f"{prefix}.tree_{idx}{suffix}*_mut_copy_numbers.tsv",
            f"{prefix}.tree_{idx}*_mut_copy_numbers.tsv",
        ])
        tot_copy_hits = collect_hits([
            f"{prefix}.tree_{idx}{suffix}_tot_copy_numbers.tsv",
            f"{prefix}.tree_{idx}{suffix}*_tot_copy_numbers.tsv",
            f"{prefix}.tree_{idx}*_tot_copy_numbers.tsv",
        ])

    if not assign_hits:
        assign_hits = collect_hits([f"*tree_{idx}_assignments.tsv"])
    if not anno_hits:
        anno_hits = collect_hits([f"*tree_{idx}.hg38_multianno.txt"])
    if not anno_hits and prefix:
        anno_hits = collect_hits([f"{prefix}.tree_*.hg38_multianno.txt"])
    if not anno_hits and prefix:
        anno_hits = sorted(glob.glob(
            os.path.join(detopt_dir, "**", f"{prefix}*.hg38_multianno.txt"),
            recursive=True,
        ))
    if not mut_copy_hits:
        mut_copy_hits = collect_hits([f"*tree_{idx}_mut_copy_numbers.tsv"])
    if not tot_copy_hits:
        tot_copy_hits = collect_hits([f"*tree_{idx}_tot_copy_numbers.tsv"])

    if not anno_hits:
        any_anno = sorted(glob.glob(os.path.join(detopt_dir, "*.hg38_multianno.txt")))
        if len(any_anno) == 1:
            anno_hits = any_anno

    assign_path = assign_hits[0] if assign_hits else None
    anno_path = anno_hits[0] if anno_hits else None
    mut_copy_path = mut_copy_hits[0] if mut_copy_hits else None
    tot_copy_path = tot_copy_hits[0] if tot_copy_hits else None
    return assign_path, anno_path, mut_copy_path, tot_copy_path


def infer_migration_counterparts(edge_path):
    """
    Given /path/foo.refined_edges_sol1[_suffix].tsv, infer:
      - nodes file: /path/foo.refined_nodes_sol1[_suffix].tsv
      - tree file:  /path/foo.tree_1[_suffix].txt
    """
    base = os.path.basename(edge_path)
    m = re.match(r"(?P<prefix>.+)\.refined_edges_sol(?P<idx>\d+)(?P<suffix>.*)\.tsv$", base)
    if not m:
        legacy = re.match(r"(?P<prefix>.+)\.refined_edges\.tsv$", base)
        if not legacy:
            return None, None
        d = os.path.dirname(edge_path)
        prefix = legacy.group("prefix")
        nodes = os.path.join(d, f"{prefix}.refined_nodes.tsv")
        tree = os.path.join(d, f"{prefix}.tree_0.txt")
        return nodes, tree
    d = os.path.dirname(edge_path)
    prefix, idx, suffix = m.group("prefix"), m.group("idx"), m.group("suffix")
    nodes = os.path.join(d, f"{prefix}.refined_nodes_sol{idx}{suffix}.tsv")
    tree = os.path.join(d, f"{prefix}.tree_{idx}{suffix}.txt")
    return nodes, tree


def looks_like_migration_edges_file(path):
    if not os.path.isfile(path) or not path.lower().endswith(".tsv"):
        return False
    try:
        df = pd.read_csv(path, sep="\t", nrows=1)
    except Exception:
        return False
    required = {"from_node", "from_site", "to_node", "to_site", "type"}
    return required.issubset(df.columns)


def build_migration_graph(edges, node_sites=None, node_to_cluster=None, site_color_map=None,
                    separate_nodes=False, node_stats=None, site_index_map=None):
    dot = graphviz.Digraph(format="pdf")
    dot.attr(rankdir="TB", nodesep="0.35", ranksep="0.45", fontsize="18")

    node_sites = node_sites or {}
    node_to_cluster = node_to_cluster or {}
    node_stats = node_stats or {}
    site_index_map = site_index_map or {}

    def make_key(node_id, site):
        base = str(node_to_cluster.get(node_id, node_id))
        site_str = "" if site is None else str(site).strip()
        if separate_nodes:
            return f"{base}|{site_str or 'NA'}"
        return base

    node_meta = {}  # key -> {node_id, site, base}
    site_refs = defaultdict(set)
    keyed_edges = []
    for e in edges:
        fn = str(e.get("from_node"))
        tn = str(e.get("to_node"))
        fs = e.get("from_site")
        ts = e.get("to_site")
        fk = make_key(fn, fs)
        tk = make_key(tn, ts)
        keyed_edges.append((fk, tk, e))
        if fs:
            site_refs[fk].add(str(fs).strip())
        if ts:
            site_refs[tk].add(str(ts).strip())
        for k, nid, st in [(fk, fn, fs), (tk, tn, ts)]:
            if k not in node_meta:
                node_meta[k] = {"node_id": nid, "site": "" if st is None else str(st).strip(),
                                 "base": str(node_to_cluster.get(nid, nid))}
            else:
                # preserve a concrete site if we previously had none
                if not node_meta[k]["site"] and st:
                    node_meta[k]["site"] = str(st).strip()

    def primary_site(key):
        info = node_meta.get(key, {})
        key_site = info.get("site")
        entries = node_sites.get(info.get("node_id")) or []
        if separate_nodes:
            if key_site:
                return key_site
            # if no explicit site, fall through to best entry
        best_site, best_signal = None, -1.0
        for ent in entries:
            presence_signal = ent.get("presence_signal")
            site = ent.get("site")
            if presence_signal is None or site in (None, ""):
                continue
            if presence_signal > best_signal:
                best_signal = presence_signal
                best_site = site
        if best_site:
            return best_site
        if site_refs.get(key):
            return sorted(site_refs[key])[0]
        return key_site or None

    def node_label(key):
        info = node_meta.get(key, {})
        node_id = info.get("node_id")
        cluster = node_to_cluster.get(node_id)
        header = f"Cluster {cluster}" if cluster else f"Node {node_id}"
        # Always show the original node id when cluster info exists
        if cluster:
            header = f"Cluster {cluster} (node {node_id})"

        key_site = info.get("site") or ""
        entries = node_sites.get(node_id) or []
        if separate_nodes and key_site:
            filtered = [ent for ent in entries if (ent.get("site") or "") == key_site]
            if filtered:
                entries = filtered
        if not entries and not site_refs.get(key):
            if key_site:
                return f"{header}\nSite: {key_site}"
            return header
        lines = [header, "Sites:"]
        if entries:
            for ent in sorted(entries, key=lambda x: x.get("site", "")):
                site = ent.get("site") or "?"
                parts = [site]
                presence_signal = ent.get("presence_signal")
                if presence_signal is not None:
                    parts.append(f"presence_signal={presence_signal:.4f}")
                obs = ent.get("observed")
                if obs not in (None, 0, 0.0):
                    parts.append(f"obs={obs:g}")
                if ent.get("is_origin"):
                    parts.append("origin")
                lines.append(" ".join(parts))
        else:
            for site in sorted(site_refs[key]):
                lines.append(site)

        stats = node_stats.get(node_id)
        if stats:
            def fmt(vals):
                if vals is None:
                    return "--"
                return ", ".join(f"{v:.2f}" for v in vals)

            idx = None
            if key_site and key_site in site_index_map:
                idx = site_index_map[key_site]

            def pick(vals):
                if vals is None or idx is None:
                    return None
                return vals[idx] if idx < len(vals) else None

            if idx is not None:
                n_val = pick(stats.get("base"))
                l_val = pick(stats.get("lineage"))
                o_val = pick(stats.get("observed"))
                lines.append(f"N: {n_val:.2f}" if n_val is not None else "N: --")
                lines.append(f"L: {l_val:.2f}" if l_val is not None else "L: --")
                if o_val is not None:
                    lines.append(f"O: {o_val:.2f}")
            else:
                lines.append(f"N: {fmt(stats.get('base'))}")
                lines.append(f"L: {fmt(stats.get('lineage'))}")
                obs = stats.get("observed")
                if obs is not None:
                    lines.append(f"O: {fmt(obs)}")

        return "\n".join(lines)

    def sort_key(val):
        base = val.split("|")[0]
        return (0, int(base)) if base.isdigit() else (1, base)

    for key in sorted(node_meta, key=sort_key):
        fill = "lightyellow"
        ps = primary_site(key)
        if ps and site_color_map and ps in site_color_map:
            fill = site_color_map[ps]
        dot.node(str(key), label=node_label(key), shape="box", style="rounded,filled", fillcolor=fill)

    for fk, tk, e in keyed_edges:
        etype = (e.get("type") or "").lower()
        color = "firebrick" if etype == "migration" else "black"
        style = "dashed" if etype == "migration" else "solid"
        parts = []
        fs, ts = e.get("from_site"), e.get("to_site")
        if fs or ts:
            parts.append(f"{fs or ''}→{ts or ''}")
        if etype:
            parts.append(etype)
        label = " ".join(parts).strip()
        dot.edge(str(fk), str(tk), color=color, fontcolor=color, style=style, label=label)

    return dot


def build_detopt_graph(data, cluster_genes, cluster_mut_counts=None, site_color_map=None,
                       gene_roles=None, include_gene_other=False):
    """
    Build a tree with cluster labels and gene lists from detopt assignments.
    cluster_genes: dict cluster -> {gene: {snv: bool, cna_states: {suffix: state}}}
    """
    dot = graphviz.Digraph(format="pdf")
    dot.attr(rankdir="TB", nodesep="0.35", ranksep="0.45", fontsize="18", ordering="out")

    node_to_cluster = {node: cl for cl, node in data["cluster_to_node"].items()}

    # nodes to render
    nodes = set()
    for p, c in data["parent_child"]:
        nodes.add(p); nodes.add(c)
    nodes.update(node_to_cluster.keys())

    def gene_lines(cl):
        genes = cluster_genes.get(str(cl), {})
        lines = []
        grouped = group_genes_by_role(genes.keys(), gene_roles, include_other=include_gene_other)
        if grouped:
            for role in ROLE_ORDER:
                if role not in grouped:
                    continue
                color = ROLE_COLORS.get(role, "black")
                lines.append(f"<FONT COLOR='{color}'><B>{role.upper()}</B></FONT>")
                for g in sorted(grouped[role]):
                    info = genes[g]
                    text = f"<FONT COLOR='{color}'>{g}</FONT>"
                    states = [val for _, val in sorted(info.get("cna_states", {}).items()) if val]
                    state_text = format_states(states)
                    if state_text:
                        text += f" {state_text}"
                    lines.append(text)
        else:
            for g in sorted(genes):
                info = genes[g]
                text = g
                states = [val for _, val in sorted(info.get("cna_states", {}).items()) if val]
                state_text = format_states(states)
                if state_text:
                    text += f" {state_text}"
                lines.append(text)
        return lines

    for n in sorted_nodes_by_cluster(nodes, node_to_cluster):
        cl = node_to_cluster.get(n)
        if str(n) == "0":
            header = "<B>Root</B>"
        else:
            header = f"<B>Cluster {cl}</B>" if cl is not None else f"<B>Node {n}</B>"

        mut_ct = cluster_mut_counts.get(str(cl), 0) if cluster_mut_counts else 0

        label_rows = [header, f"Mutations: {mut_ct}"]
        if cl is not None:
            gene_list = gene_lines(cl)
            if gene_list:
                label_rows.append("__DIV__")
                label_rows.extend(gene_list)

        # left-aligned table with a small left margin spacer
        row_parts = []
        for row in label_rows:
            if row == "__DIV__":
                row_parts.append("<TR><TD ALIGN='LEFT'><FONT COLOR='gray'>────────</FONT></TD></TR>")
            else:
                row_parts.append(f"<TR><TD ALIGN='LEFT'>{row}</TD></TR>")
        rows_html = "".join(row_parts)
        label_html = (
            "<<TABLE BORDER='0' CELLBORDER='0' CELLSPACING='0' CELLPADDING='4' ALIGN='RIGHT'>"
            "<TR><TD WIDTH='8'></TD></TR>"
            f"{rows_html}"
            "</TABLE>>"
        )
        fill = "lightyellow"
        dot.node(str(n), label=label_html, shape="box", style="rounded,filled", fillcolor=fill)

    for p, c in sorted_parent_child_edges(data["parent_child"], node_to_cluster):
        dot.edge(p, c)

    return dot


def build_detopt_graph_main_figure(data, cluster_genes, cluster_mut_counts=None,
                                   gene_roles=None, include_gene_other=False):
    """
    Alternative detopt-style plot roughly matching DETOPT demo (circular nodes with side gene labels).
    """
    dot = graphviz.Digraph(format="pdf")
    dot.attr(rankdir="TB", nodesep="0.45", ranksep="0.6", fontsize="18", ordering="out")

    node_to_cluster = {node: cl for cl, node in data["cluster_to_node"].items()}

    def gene_rows(cl):
        genes = cluster_genes.get(str(cl), {})
        rows = []
        grouped = group_genes_by_role(genes.keys(), gene_roles, include_other=include_gene_other)
        if grouped:
            for role in ROLE_ORDER:
                if role not in grouped:
                    continue
                color = ROLE_COLORS.get(role, "black")
                rows.append(
                    f"<TR><TD ALIGN='LEFT'><FONT COLOR='{color}'><B>{role.upper()}</B></FONT></TD></TR>"
                )
                for g in sorted(grouped[role]):
                    info = genes[g]
                    text = f"<FONT COLOR='{color}'><I>{g}</I></FONT>"
                    states = [val for _, val in sorted(info.get("cna_states", {}).items()) if val]
                    state_text = format_states(states)
                    if state_text:
                        text += f" <FONT COLOR='purple'>{state_text}</FONT>"
                    rows.append(f"<TR><TD ALIGN='LEFT'>{text}</TD></TR>")
        else:
            for g in sorted(genes):
                info = genes[g]
                text = f"<I>{g}</I>"
                states = [val for _, val in sorted(info.get("cna_states", {}).items()) if val]
                state_text = format_states(states)
                if state_text:
                    text += f" <FONT COLOR='purple'>{state_text}</FONT>"
                rows.append(f"<TR><TD ALIGN='LEFT'>{text}</TD></TR>")
        return "".join(rows)

    nodes = set()
    for p, c in data["parent_child"]:
        nodes.add(p); nodes.add(c)
    nodes.update(node_to_cluster.keys())

    for n in sorted_nodes_by_cluster(nodes, node_to_cluster):
        cl = node_to_cluster.get(n)
        mut_ct = cluster_mut_counts.get(str(cl), 0) if cluster_mut_counts else 0
        header = "ROOT" if str(n) == "0" else (str(cl) if cl is not None else str(n))
        fill = "dodgerblue" if str(n) == "0" else "white"
        fontcolor = "white" if str(n) == "0" else "black"

        gene_html = gene_rows(cl) if cl is not None else ""
        xlabel = ""
        if gene_html:
            xlabel = (
                "<"
                "<TABLE BORDER='0' CELLBORDER='0' CELLSPACING='0' CELLPADDING='1'>"
                f"{gene_html}"
                "</TABLE>"
                ">"
            )

        dot.node(
            str(n),
            label=f"{header}\\nMutations: {mut_ct}",
            shape="circle",
            style="filled",
            fillcolor=fill,
            fontcolor=fontcolor,
            xlabel=xlabel,
            fixedsize="false",
            width="0.85",
            height="0.85",
            penwidth="1.5",
        )

    for p, c in sorted_parent_child_edges(data["parent_child"], node_to_cluster):
        dot.edge(p, c, arrowsize="0.8")

    return dot


# ---------------- graph build ----------------

def _float_equal(a_str, b_str, tol=1e-12):
    try:
        return abs(float(a_str) - float(b_str)) <= tol
    except Exception:
        return False


def build_graph(data, gene_info=None, node_to_site=None, site_color_map=None, sample_site_map=None,
               gene_wrap=GENE_LINE_WRAP, include_gene_other=False, transitions=None, allowed_genes=None,
               with_node_info=False, show_presence=False, presence_by_node=None,
               presence_site_order=None, presence_origin_sites=None, extra_title_lines=None):
    dot = graphviz.Digraph(format="pdf")
    dot.attr(rankdir="TB", nodesep="0.35", ranksep="0.45", fontsize="20", ordering="out")

    # children
    children = defaultdict(list)
    for p, c in data["parent_child"]:
        children[p].append(c)

    # mutation counts
    mut_count = {}
    for mut, n in data["mutation_to_node"].items():
        mut_count[n] = mut_count.get(n, 0) + 1
    total_mut = len(data["mutation_to_node"])

    # average observed per node
    obs_sum, obs_cnt = {}, {}
    for mut, (obs, _, _) in data["mutation_info"].items():
        n = data["mutation_to_node"].get(mut)
        if not n:
            continue
        try:
            vals = [float(x.strip()) for x in obs.split(",") if x.strip()]
        except Exception:
            continue
        if n not in obs_sum:
            obs_sum[n] = [0.0] * len(vals)
            obs_cnt[n] = 0
        elif len(obs_sum[n]) != len(vals):
            continue
        for i, v in enumerate(vals):
            obs_sum[n][i] += v
        obs_cnt[n] += 1

    nodes = set()
    for p, c in data["parent_child"]:
        nodes.add(p); nodes.add(c)
    nodes.update(data["node_freqs"].keys())
    nodes.update(data["cluster_to_node"].values())
    nodes.update(obs_sum.keys())

    avg_obs = {n: (tuple(v / obs_cnt[n] for v in obs_sum[n]) if n in obs_cnt and obs_cnt[n] else None) for n in nodes}
    node_to_cluster = {node: cl for cl, node in data["cluster_to_node"].items()}

    # lineage totals
    total_freqs = {}
    def acc(n):
        if n in total_freqs:
            return total_freqs[n]
        base = data["node_freqs"].get(n)
        total = list(base) if base else []
        for ch in children.get(n, []):
            total = add_freq_lists(total, acc(ch))
        total_freqs[n] = total
        return total
    for n in sorted_nodes_by_cluster(nodes, node_to_cluster):
        acc(n)

    # samples width
    sample_count = 0
    if total_freqs:
        sample_count = max((len(v) for v in total_freqs.values() if v), default=0)
    if sample_count == 0:
        sample_count = max((len(v) for v in data["node_freqs"].values() if v), default=0)
    sample_items = sample_site_items_for_width(sample_site_map, sample_count)
    presence_site_order = list(presence_site_order or [])
    presence_by_node = presence_by_node or {}
    presence_origin_sites = presence_origin_sites or node_to_site or {}

    def presence_bits_for_node(node_id):
        site_map = presence_by_node.get(str(node_id))
        if not site_map:
            return None

        sample_sites = [site for _, site in sample_items if site]
        if sample_sites:
            ordered_sites = [site for _, site in sample_items]
        else:
            ordered_sites = list(presence_site_order)

        if not ordered_sites:
            return None
        return effective_presence_bits(node_id, site_map, ordered_sites, presence_origin_sites)

    # ---------- top label: centered title, legend top-right ----------
    tree_score = data["tree_score"] or "N/A"
    raw_score = data.get("raw_tree_score") or "N/A"
    show_raw = (raw_score != "N/A") and (not _float_equal(tree_score, raw_score))

    def centered_stats():
        parts = [f"<B>Tree Score:</B> {tree_score}"]
        if show_raw:
            parts.append(f"<B>Tree Score (no reward):</B> {raw_score}")
        parts.extend(extra_title_lines or [])
        parts.append(f"<B>Total mutations:</B> {total_mut}")
        # One cell, multiple stacked rows centered
        inner = "".join(f"<TR><TD ALIGN='CENTER'>{line}</TD></TR>" for line in parts)
        return f"<TABLE BORDER='0' CELLBORDER='0' CELLSPACING='0' CELLPADDING='2'>{inner}</TABLE>"

    def legend_table():
        if not sample_site_map:
            return ""
        items = ordered_sample_site_items(sample_site_map)
        header = (
            "<TR>"
            "<TD ALIGN='CENTER'><B>#</B></TD>"
            "<TD ALIGN='LEFT'><B>Sample</B></TD>"
            "<TD ALIGN='LEFT'><B>Site</B></TD>"
            "<TD ALIGN='CENTER'><B>Color</B></TD>"
            "</TR>"
        )
        body = []
        for i, (sample, site) in enumerate(items, start=1):
            site_str = site if site else ""
            color = site_color_map.get(site_str, "white") if site_color_map else "white"
            body.append(
                "<TR>"
                f"<TD ALIGN='CENTER'>{i}</TD>"
                f"<TD ALIGN='LEFT'>{sample}</TD>"
                f"<TD ALIGN='LEFT' BGCOLOR='{color}'>{site_str}</TD>"
                f"<TD WIDTH='18' HEIGHT='12' BGCOLOR='{color}'></TD>"
                "</TR>"
            )
        return ("<TABLE BORDER='1' CELLBORDER='1' CELLSPACING='0' CELLPADDING='6'>"
                "<TR><TD COLSPAN='4' ALIGN='LEFT'><B>Sample → Site</B></TD></TR>"
                f"{header}{''.join(body)}"
                "</TABLE>")

    stats = centered_stats()
    legend = legend_table()

    if legend:
        # Row 1: centered stats across both columns
        # Row 2: big left spacer (fixed width) + right-aligned legend
        # The spacer width pushes the legend to the far right of the label.
        label = (
            "<<TABLE BORDER='0' CELLBORDER='0' CELLSPACING='0' WIDTH='100%'>"
            "<TR><TD ALIGN='CENTER' VALIGN='TOP' COLSPAN='2'>"
            f"{stats}"
            "</TD></TR>"
            "<TR>"
            "<TD WIDTH='900'></TD>"  # <-- spacer; bump this if you want the legend even farther right
            f"<TD ALIGN='RIGHT' VALIGN='TOP'>{legend}</TD>"
            "</TR>"
            "</TABLE>>"
        )
    else:
        label = (
            "<<TABLE BORDER='0' CELLBORDER='0' CELLSPACING='0' WIDTH='100%'>"
            "<TR><TD ALIGN='CENTER' VALIGN='TOP'>"
            f"{stats}"
            "</TD></TR>"
            "</TABLE>>"
        )

    dot.attr(labelloc="t", label=label)

    # ---------- nodes ----------
    for n in sorted_nodes_by_cluster(nodes, node_to_cluster):
        site = node_to_site.get(str(n)) if node_to_site else None
        site = site or None
        fill = "lightyellow"
        if site and site_color_map and site in site_color_map:
            fill = site_color_map[site]

        if n == "0":
            dot.node("0", label="ROOT", shape="box", style="rounded,filled", fillcolor=fill)
            continue

        base = data["node_freqs"].get(n)
        lin = total_freqs.get(n)
        obs = avg_obs.get(n)

        widest = max([len(v) for v in (base or [], lin or [], obs or [])] + [sample_count or 1])

        def hrow(txt): return f"<TR><TD COLSPAN='{widest + 1}' ALIGN='LEFT'><B>{txt}</B></TD></TR>"

        def frow(name, arr):
            cells = [f"<TD ALIGN='LEFT'><B>{name}</B></TD>"]
            vals = list(arr) if arr else []
            for i in range(widest):
                cells.append(f"<TD ALIGN='RIGHT'>{vals[i]:.2f}</TD>" if i < len(vals) else "<TD ALIGN='RIGHT'>--</TD>")
            return f"<TR>{''.join(cells)}</TR>"

        def brow(name, arr):
            cells = [f"<TD ALIGN='LEFT'><B>{name}</B></TD>"]
            vals = list(arr) if arr else []
            for i in range(widest):
                if i < len(vals):
                    try:
                        bit = 1 if float(vals[i]) > 0 else 0
                    except Exception:
                        bit = "--"
                    cells.append(f"<TD ALIGN='RIGHT'>{bit}</TD>")
                else:
                    cells.append("<TD ALIGN='RIGHT'>--</TD>")
            return f"<TR>{''.join(cells)}</TR>"

        rows = []
        if n in node_to_cluster:
            rows.append(hrow(f"Cluster: {node_to_cluster[n]}"))
        if with_node_info:
            rows.append(hrow(f"Node: {n}"))
        if site:
            rows.append(hrow(f"Site: {site}"))
        rows.append(hrow(f"Mutations: {mut_count.get(n, 0)}"))
        if show_presence:
            presence_bits = presence_bits_for_node(n)
            if presence_bits is not None:
                rows.append(brow("Node &gt;0", presence_bits))
            elif sample_items:
                rows.append(brow("Node &gt;0", base))

        if base or lin or obs:
            rows.append(frow("Inf Node", base))
            rows.append(frow("Inf Lin", lin))
            if obs is not None:
                rows.append(frow("Obs Lin", obs))
        else:
            rows.append(hrow("No frequency data"))

        # gene note (inline under node content)
        if n in node_to_cluster and gene_info:
            mutated = set()
            for mut, nn in data["mutation_to_node"].items():
                if nn == n:
                    chrom, pos = parse_mutation_coordinate(mut)
                    if chrom is None:
                        continue
                    for gene, info in gene_info.items():
                        if allowed_genes is not None and gene.lower() not in allowed_genes:
                            continue
                        if info["chrom"] == chrom and info["start"] <= pos <= info["end"]:
                            mutated.add(gene)
            if mutated:
                role_to_genes = {}
                for g in mutated:
                    roles = classify_gene_roles(gene_info[g].get("role", ""), include_other=include_gene_other)
                    for r in roles:
                        role_to_genes.setdefault(r, set()).add(g)
                if role_to_genes:
                    # separator between frequency rows and gene annotations
                    rows.append(
                        f"<TR><TD COLSPAN='{widest + 1}' HEIGHT='1' BGCOLOR='black'></TD></TR>"
                    )
                    ordered_roles = [r for r in ROLE_ORDER if r in role_to_genes]
                    for role in ordered_roles:
                        color = ROLE_COLORS.get(role, "black")
                        gene_lines = wrap_role_gene_lines(role, role_to_genes[role], color, max_chars=gene_wrap)
                        if not gene_lines:
                            continue
                        gene_span = max(1, widest)
                        # first line includes the role label in its own cell, continuation lines leave it blank
                        rows.append(
                            "<TR>"
                            f"<TD ALIGN='LEFT'><FONT COLOR='{color}'><B>{role.upper()}</B></FONT></TD>"
                            f"<TD ALIGN='LEFT' COLSPAN='{gene_span}'><FONT COLOR='{color}'>{gene_lines[0]}</FONT></TD>"
                            "</TR>"
                        )
                        for cont in gene_lines[1:]:
                            rows.append(
                                "<TR>"
                                "<TD></TD>"
                                f"<TD ALIGN='LEFT' COLSPAN='{gene_span}'><FONT COLOR='{color}'>{cont}</FONT></TD>"
                                "</TR>"
                            )

        dot.node(n,
                 label="<<TABLE BORDER='0' CELLBORDER='0' CELLSPACING='4'>" + "".join(rows) + "</TABLE>>",
                 shape="box", style="rounded,filled", fillcolor=fill)

    # edges
    for p, c in sorted_parent_child_edges(data["parent_child"], node_to_cluster):
        dot.edge(p, c)

    # optional site transitions as self-loops
    if transitions:
        for n, tr_list in transitions.items():
            grouped = defaultdict(int)
            has_migration = False
            for tr in tr_list:
                etype = (tr.get("type") or "").lower()
                if etype == "migration":
                    has_migration = True
                fs, ts = tr.get("from_site"), tr.get("to_site")
                grouped[(fs, ts, etype)] += 1

            lines = []
            def sort_key(item):
                fs, ts, etype = item[0]
                return (fs or "", ts or "", etype or "")

            for (fs, ts, etype), count in sorted(grouped.items(), key=sort_key):
                arrow = f"{fs or ''}→{ts or ''}"
                line = f"{arrow} {etype}".strip()
                if count > 1:
                    line += f" (x{count})"
                lines.append(line)

            label = "\\n".join(lines)
            color = "firebrick" if has_migration else "black"
            style = "dashed" if has_migration else "dotted"
            dot.edge(str(n), str(n), color=color, fontcolor=color, style=style,
                     label=label, labeldistance="2.2", penwidth="2")

    return dot


def render_migration_plots(edge_files, args, allowed_genes=None):
    for edge_file in edge_files:
        edges = parse_migration_edges(edge_file)
        if not edges:
            continue

        inferred_nodes, inferred_tree = infer_migration_counterparts(edge_file)
        nodes_path = args.refined_nodes or inferred_nodes
        tree_path = args.refined_tree or inferred_tree

        # Fallback: if we failed to infer, try matching by prefix before ".refined_edges"
        if (not tree_path) or (not os.path.exists(tree_path)):
            base = os.path.basename(edge_file)
            prefix = base.split(".refined_edges")[0]
            candidate_trees = sorted(glob.glob(os.path.join(os.path.dirname(edge_file), f"{prefix}.tree_*.txt")))
            if candidate_trees:
                tree_path = candidate_trees[0]

        site_index_map = {}
        if nodes_path and not os.path.exists(nodes_path):
            sys.stderr.write(f"Warning: refined nodes file not found: {nodes_path}\n")
        node_sites = parse_migration_nodes(nodes_path) if nodes_path and os.path.exists(nodes_path) else {}
        if nodes_path and os.path.exists(nodes_path):
            order = site_order_from_migration_nodes(nodes_path)
            site_index_map = {s: i for i, s in enumerate(order)}

        node_to_cluster = {}
        node_stats = {}
        if tree_path:
            if os.path.exists(tree_path):
                tdata = parse_tree_file(tree_path)
                node_to_cluster = {str(node): str(cluster) for cluster, node in tdata["cluster_to_node"].items()}
                node_stats = compute_node_stats(tdata)
            else:
                sys.stderr.write(f"Warning: cluster-to-node map not found: {tree_path}\n")

        edge_sites = {s for e in edges for s in (e.get("from_site"), e.get("to_site")) if s}
        node_site_names = {ent["site"] for ent_list in node_sites.values() for ent in ent_list if ent.get("site")}
        node_site_order = site_order_from_migration_nodes(nodes_path) if nodes_path and os.path.exists(nodes_path) else []
        site_color_map = make_site_color_map(
            edge_sites,
            node_site_names,
            site_order=node_site_order,
            override_colors=getattr(args, "site_color_overrides", None),
        ) if (edge_sites or node_site_names or node_site_order) else None

        dot = build_migration_graph(edges, node_sites=node_sites, node_to_cluster=node_to_cluster,
                               site_color_map=site_color_map, separate_nodes=args.separate_nodes,
                               node_stats=node_stats, site_index_map=site_index_map)

        if args.output:
            base, ext = os.path.splitext(args.output)
            if ext and ext.lower() != f".{dot.format}":
                sys.stderr.write(f"Warning: ignoring output extension '{ext}' in favor of '.{dot.format}'.\n")
            out_no_ext = base if ext else args.output
        else:
            base = os.path.splitext(os.path.basename(edge_file))[0]
            out_no_ext = os.path.join(args.output_dir or os.path.dirname(edge_file), base)

        out_path = dot.render(filename=out_no_ext, cleanup=True)
        msg = f"Saved migration plot to {out_path}"
        if nodes_path and os.path.exists(nodes_path):
            msg += f"  (nodes: {os.path.basename(nodes_path)})"
        if tree_path and os.path.exists(tree_path):
            msg += f"  (clusters: {os.path.basename(tree_path)})"
        print(msg)


# --------------- main ---------------

def main():
    parser = argparse.ArgumentParser(
        description=("Plot CITUPm tree(s) with node coloring by site. Legend sits at the top-right; "
                     "title/scores centered; 'no reward' hidden if identical; no 'Samples:' line.")
    )
    parser.add_argument("tree_path", help="*.tree_*.txt file or a directory containing such files")
    parser.add_argument("--gene-info", help="Optional gene annotation TSV")
    parser.add_argument("--gene-wrap", type=int, default=GENE_LINE_WRAP,
                        help="Approx max characters per gene line (smaller makes nodes narrower)")
    parser.add_argument("--include-gene-other", action="store_true",
                        help="Show genes without 'Role in Cancer' as OTHER")
    parser.add_argument("--sample-sites", help="Optional TSV mapping: sample\\tsite (adds colored, numbered legend)")
    parser.add_argument("--site-colors",
                        help="Optional TSV mapping: site\\tcolor. Use #RRGGBB hex colors or Graphviz color names.")
    parser.add_argument("--output", help="Destination PDF path (single-file mode only)")
    parser.add_argument("--output-dir", help="Directory for generated PDFs when plotting multiple files")
    parser.add_argument("--refined-nodes", help="Optional refined_nodes_sol*.tsv to annotate migration plots")
    parser.add_argument("--refined-tree", help="Optional tree_* file to map node ids to cluster ids for migration plots")
    parser.add_argument("--separate-nodes", action="store_true",
                        help="Treat each node+site combination as a distinct node in migration plots")
    parser.add_argument("--with-site-transitions", action="store_true",
                        help="Overlay site transitions (from matching migration_edges files) as self-loops on tree plots")
    parser.add_argument("--migration-dag-coloring-index", type=int,
                        help="Use embedded DAG_COLORING_INDEX <k> from the tree file for node site coloring")
    parser.add_argument("--all-migration-dag-colorings", action="store_true",
                        help="Render one tree plot per embedded DAG_COLORING_INDEX entry")
    parser.add_argument("--all-solutions-all-dags", action="store_true",
                        help="Render every primary tree solution and every embedded CITUPm DAG coloring, with site transitions")
    parser.add_argument("--include-best-migration-copy", action="store_true",
                        help="Also plot *_best_migration tree copies when plotting a directory")
    parser.add_argument("--with_node_info", "--with-node-info", dest="with_node_info", action="store_true",
                        help='Add "Node: <id>" in each non-root node label on main tree plots')
    parser.add_argument("--detopt", help="Directory containing detopt assignments/multianno files (adds extra gene-only tree plot)")
    parser.add_argument("--intersection", action="store_true",
                        help="Only show genes present in --gene-info (symbol or any synonym)")
    args = parser.parse_args()

    if not os.path.exists(args.tree_path):
        parser.error(f"Path '{args.tree_path}' does not exist.")

    if args.output_dir:
        os.makedirs(args.output_dir, exist_ok=True)
    if args.all_solutions_all_dags:
        args.all_migration_dag_colorings = True
        args.with_site_transitions = True
    args.site_color_overrides = parse_site_colors(args.site_colors) if args.site_colors else {}

    # Decide whether we are plotting refined migration TSVs or standard tree files.
    migration_edge_files = []
    tree_files = []
    if os.path.isdir(args.tree_path):
        tree_files = discover_tree_files(args.tree_path, include_best_migration=args.include_best_migration_copy)
        if not tree_files:
            migration_edge_files = sorted(glob.glob(os.path.join(args.tree_path, "*.refined_edges*.tsv")))
    else:
        if looks_like_migration_edges_file(args.tree_path):
            migration_edge_files = [args.tree_path]
        else:
            tree_files = [args.tree_path]

    if migration_edge_files:
        if args.output and len(migration_edge_files) > 1:
            parser.error("--output is only valid when plotting a single file.")
        render_migration_plots(migration_edge_files, args)
        return

    if not tree_files:
        parser.error(f"No '*.tree_*.txt' files found under '{args.tree_path}'.")

    if args.output and len(tree_files) > 1:
        parser.error("--output is only valid when plotting a single tree file.")
    if args.migration_dag_coloring_index is not None and args.all_migration_dag_colorings:
        parser.error("--migration-dag-coloring-index and --all-migration-dag-colorings cannot be used together.")

    gene_info = parse_gene_info(args.gene_info, include_no_role=args.include_gene_other) if args.gene_info else None
    gene_roles = build_gene_role_map(gene_info) if gene_info else None
    allowed_genes = None
    if args.intersection:
        if gene_info:
            allowed_genes = allowed_gene_set(gene_info)
        else:
            sys.stderr.write("--intersection requested but --gene-info not provided; ignoring intersection filter.\n")
    explicit_sample_site_map = parse_sample_sites(args.sample_sites) if args.sample_sites else None

    tree_data_cache = {}
    global_site_order = []
    seen_global_sites = set()

    def add_global_sites(sites):
        for site in sites or []:
            if site and site not in seen_global_sites:
                seen_global_sites.add(site)
                global_site_order.append(site)

    for tree_file in tree_files:
        data = parse_tree_file(tree_file)
        tree_data_cache[tree_file] = data
        add_global_sites(data.get("m_site_order", []))
        add_global_sites(data.get("origin_sites", {}).values())
        for site_map in (data.get("dag_colorings") or {}).values():
            add_global_sites(site_map.values())
        for site_map in (data.get("m_presence") or {}).values():
            add_global_sites(site_map.keys())
        for edge in data.get("migration_edges", []):
            add_global_sites([edge.get("from_site"), edge.get("to_site")])

    if explicit_sample_site_map:
        add_global_sites(site for _, site in ordered_sample_site_items(explicit_sample_site_map))

    for tree_file in tree_files:
        data = tree_data_cache.get(tree_file) or parse_tree_file(tree_file)
        sample_site_map = explicit_sample_site_map
        if sample_site_map is None:
            inferred_sample_sites = find_sample_sites_file_for_tree(tree_file)
            if inferred_sample_sites:
                sample_site_map = parse_sample_sites(inferred_sample_sites)
        dag_colorings = data.get("dag_colorings") or {}
        if args.all_migration_dag_colorings:
            selected_colorings = sorted(dag_colorings.keys())
            if not selected_colorings:
                sys.stderr.write(f"Warning: no embedded CITUPm DAG colorings found for {tree_file}\n")
                selected_colorings = [None]
        else:
            selected_colorings = [args.migration_dag_coloring_index]

        for requested_coloring in selected_colorings:
            node_to_site, site_color_map = None, None
            site_file = find_site_labels_file_for_tree(tree_file)
            selected_dag_coloring = None
            dag_title_lines = []
            header_values = data.get("header_values") or {}
            if header_values.get("REFINED_SOLUTION_INDEX") is not None:
                m_title = f"<B>Refined solution:</B> {header_values.get('REFINED_SOLUTION_INDEX')}"
                if header_values.get("REFINED_SOLUTION_COUNT") is not None:
                    m_title += f" / {header_values.get('REFINED_SOLUTION_COUNT')}"
                dag_title_lines.append(m_title)
            if requested_coloring is not None:
                node_to_site, selected_dag_coloring = choose_origin_sites(data, requested_coloring)
                if node_to_site is None:
                    avail = sorted(dag_colorings.keys())
                    sys.stderr.write(
                        f"Warning: CITUPm DAG coloring index {requested_coloring} not found for {tree_file}"
                        + (f"; available: {avail}\n" if avail else "; no embedded DAG colorings present\n")
                    )
                    continue
                stats = data.get("dag_coloring_stats", {}).get(selected_dag_coloring, {})
                dag_title = f"<B>CITUPm DAG coloring:</B> {selected_dag_coloring}"
                dag_bits = []
                if stats.get("same_origin_score") is not None:
                    dag_bits.append(f"same-origin {stats.get('same_origin_score')}")
                if stats.get("migration_edge_count") is not None:
                    dag_bits.append(f"migration edges {stats.get('migration_edge_count')}")
                if stats.get("dag") is not None:
                    dag_bits.append(f"DAG {stats.get('dag')}")
                if dag_bits:
                    dag_title += " (" + ", ".join(dag_bits) + ")"
                dag_title_lines.append(dag_title)
            if node_to_site:
                site_file = None
            elif site_file:
                node_to_site = parse_site_labels(site_file)
            elif data.get("origin_sites"):
                # Fallback for migration-augmented tree files that already embed node->site.
                node_to_site = data["origin_sites"]

            transitions = None
            transition_sites = set()
            if args.with_site_transitions:
                m_presence = data.get("m_presence") or {}
                m_site_order = data.get("m_site_order") or []

                edges = []
                transitions = None
                if selected_dag_coloring is not None and node_to_site and m_presence:
                    transitions, transition_sites = derive_migration_transitions(node_to_site, m_presence, m_site_order)
                else:
                    edges_path = find_migration_edges_for_tree(tree_file)
                    if edges_path and os.path.exists(edges_path):
                        edges = parse_migration_edges(edges_path)
                    if not edges:
                        edges = data.get("migration_edges", [])

                if edges:
                    transitions = defaultdict(list)
                    for e in edges:
                        n_from = e.get("from_node")
                        n_to = e.get("to_node")
                        from_site = e.get("from_site")
                        to_site = e.get("to_site")
                        if n_from != n_to:
                            continue  # only self transitions requested
                        transitions[str(n_from)].append({
                            "from_site": from_site,
                            "to_site": to_site,
                            "type": e.get("type"),
                        })
                        for s in (from_site, to_site):
                            if s:
                                transition_sites.add(s)
                elif not transitions:
                    sys.stderr.write(f"Warning: no transition edges found in migration_edges or tree file for {tree_file}\n")

            sites_in_nodes = set(node_to_site.values()) if node_to_site else set()
            sites_in_samples = set(sample_site_map.values()) if sample_site_map else set()
            if sites_in_nodes or sites_in_samples or transition_sites:
                plot_site_order = list(global_site_order)
                plot_site_order.extend(data.get("m_site_order") or [])
                if sample_site_map:
                    plot_site_order.extend(site for _, site in ordered_sample_site_items(sample_site_map))
                site_color_map = make_site_color_map(
                    sites_in_nodes,
                    sites_in_samples,
                    transition_sites,
                    site_order=plot_site_order,
                    override_colors=args.site_color_overrides,
                )

            gene_wrap = max(10, args.gene_wrap)  # avoid degenerate widths
            dot = build_graph(data, gene_info, node_to_site, site_color_map, sample_site_map,
                              gene_wrap=gene_wrap, include_gene_other=args.include_gene_other,
                              transitions=transitions, allowed_genes=allowed_genes,
                              with_node_info=args.with_node_info, show_presence=args.with_site_transitions,
                              presence_by_node=data.get("m_presence"),
                              presence_site_order=data.get("m_site_order"),
                              presence_origin_sites=node_to_site,
                              extra_title_lines=dag_title_lines)

            if args.output:
                base, ext = os.path.splitext(args.output)
                if ext and ext.lower() != f".{dot.format}":
                    sys.stderr.write(f"Warning: ignoring output extension '{ext}' in favor of '.{dot.format}'.\n")
                out_no_ext = base if ext else args.output
                if selected_dag_coloring is not None:
                    out_no_ext = f"{out_no_ext}.migration_dag_coloring_{selected_dag_coloring}"
            else:
                base = os.path.splitext(os.path.basename(tree_file))[0]
                if selected_dag_coloring is not None:
                    base += f".migration_dag_coloring_{selected_dag_coloring}"
                out_no_ext = os.path.join(args.output_dir or os.path.dirname(tree_file), base)

            out_path = dot.render(filename=out_no_ext, cleanup=True)
            print(f"Saved tree plot to {out_path}" + (f"  (sites: {os.path.basename(site_file)})" if site_file else ""))

            # Optional detopt gene-only plot
            if args.detopt:
                assign_path, anno_path, mut_copy_path, tot_copy_path = find_detopt_files(args.detopt, tree_file)
                if assign_path and anno_path:
                    cluster_genes, cluster_mut_counts = build_detopt_cluster_genes(
                        assign_path,
                        anno_path,
                        allowed_genes=allowed_genes,
                        mut_copy_path=mut_copy_path,
                        tot_copy_path=tot_copy_path,
                    )
                    if cluster_genes:
                        detopt_dot = build_detopt_graph(
                            data,
                            cluster_genes,
                            cluster_mut_counts=cluster_mut_counts,
                            site_color_map=site_color_map,
                            gene_roles=gene_roles,
                            include_gene_other=args.include_gene_other,
                        )
                        detopt_out = f"{out_no_ext}_detopt"
                        detopt_path = detopt_dot.render(filename=detopt_out, cleanup=True)
                        print(f"Saved detopt tree plot to {detopt_path} (genes from {os.path.basename(assign_path)})")

                        # additional detopt-style figure
                        detopt_main = build_detopt_graph_main_figure(
                            data,
                            cluster_genes,
                            cluster_mut_counts=cluster_mut_counts,
                            gene_roles=gene_roles,
                            include_gene_other=args.include_gene_other,
                        )
                        detopt_main_out = f"{out_no_ext}_detopt_main"
                        detopt_main_path = detopt_main.render(filename=detopt_main_out, cleanup=True)
                        print(f"Saved detopt main-style plot to {detopt_main_path}")
                    else:
                        sys.stderr.write(f"No detopt genes found for {tree_file}\n")
                else:
                    sys.stderr.write(f"Detopt files not found for {tree_file} under {args.detopt}\n")


if __name__ == "__main__":
    main()
