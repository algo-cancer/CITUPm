#!/usr/bin/env python3
"""
Render CITUPm migration and comigration graphs from CITUPm sidecar files.

Inputs are the *.refined_edges_sol*.tsv files emitted by citupm.py.
For each edge file this script writes:
  - *.migration_graph.dot/pdf: expanded node-site migration graph.
  - *.comigration_graph.dot/pdf: site-level graph with direct migration counts.

For each tree topology (solN), it also writes:
  - *.comigration_conjoined_solN.dot/pdf: site-level graph summed across all
    refined trees for that topology.
  - matching *.tsv count summaries for the comigration edges.
"""

import argparse
import csv
import glob
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict


SITE_COLORS = [
    "#add8e6", "#98fb98", "#E69F00", "#e5a79e", "#a5b1eb",
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

CLONE_SITE_COLORS = [
    "#e41a1c", "#33a02c", "#377eb8", "#984ea3", "#ff7f00", "#a65628",
    "#f781bf", "#008b8b", "#666666", "#8c8b00",
]

EDGE_RE = re.compile(r"(?P<prefix>.+)\.refined_edges_sol(?P<idx>\d+)(?P<suffix>.*)\.tsv$")


def natural_key(value):
    parts = re.split(r"(\d+)", str(value))
    return [int(p) if p.isdigit() else p.lower() for p in parts]


def typed_natural_key(value):
    parts = re.split(r"(\d+)", str(value))
    return [(0, int(p)) if p.isdigit() else (1, p.lower()) for p in parts if p != ""]


def dot_quote(value):
    text = "" if value is None else str(value)
    text = text.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n")
    return f'"{text}"'


def attr_text(attrs):
    if not attrs:
        return ""
    return " [" + ", ".join(f"{k}={dot_quote(v)}" for k, v in attrs.items() if v is not None) + "]"


def html_escape(value):
    text = "" if value is None else str(value)
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def site_legend_html(site_color_map):
    if not site_color_map:
        return ""
    sites = sorted(site_color_map, key=natural_key)
    rows_per_col = 3
    columns = [sites[i:i + rows_per_col] for i in range(0, len(sites), rows_per_col)]
    lines = []
    for row_idx in range(rows_per_col):
        entries = []
        for col in columns:
            if row_idx < len(col):
                site = col[row_idx]
                color = html_escape(site_color_map.get(site, "lightyellow"))
                site_label = html_escape(site)
                entries.append(
                    f"<FONT COLOR='{color}' POINT-SIZE='16'>&#9632;</FONT>"
                    f"<FONT POINT-SIZE='12'> {site_label}</FONT>"
                )
        if entries:
            lines.append("&#160;&#160;&#160;&#160;".join(entries))
    return (
        "<TABLE BORDER='0' CELLBORDER='0' CELLSPACING='0' CELLPADDING='0'>"
        + "".join(f"<TR><TD ALIGN='LEFT'>{line}</TD></TR>" for line in lines)
        + "</TABLE>"
    )


def graph_label_html(title, site_color_map, note=None):
    title_html = html_escape(title)
    legend = site_legend_html(site_color_map)
    if not legend:
        return f"<{title_html}>"
    note_row = f"<TR><TD ALIGN='LEFT'><FONT POINT-SIZE='10'>{html_escape(note)}</FONT></TD></TR>" if note else ""
    return (
        "<<TABLE BORDER='0' CELLBORDER='0' CELLSPACING='0' CELLPADDING='0'>"
        f"<TR><TD ALIGN='LEFT'><FONT POINT-SIZE='18'>{title_html}</FONT></TD></TR>"
        f"<TR><TD ALIGN='LEFT'>{legend}</TD></TR>"
        + note_row +
        '<TR><TD HEIGHT="12"></TD></TR>' +
        "</TABLE>>"
    )


def clean_node_id(value):
    text = "" if value is None else str(value).strip()
    if not text:
        return ""
    try:
        num = float(text)
        if num.is_integer():
            return str(int(num))
    except Exception:
        pass
    return text


def parse_edges(path):
    required = {"from_node", "from_site", "to_node", "to_site", "type"}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = set(reader.fieldnames or [])
        if not required.issubset(fields):
            raise ValueError(f"{path} must contain columns: {', '.join(sorted(required))}")
        rows = []
        for row in reader:
            rows.append({
                "from_node": clean_node_id(row.get("from_node")),
                "from_site": (row.get("from_site") or "").strip(),
                "to_node": clean_node_id(row.get("to_node")),
                "to_site": (row.get("to_site") or "").strip(),
                "type": (row.get("type") or "").strip().lower(),
            })
        return rows


def parse_site_order(path):
    if not path or not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            return []
        col = "site" if "site" in reader.fieldnames else reader.fieldnames[-1]
        seen, order = set(), []
        for row in reader:
            site = (row.get(col) or "").strip()
            if site and site not in seen:
                seen.add(site)
                order.append(site)
        return order


def parse_site_colors(path):
    if not path:
        return {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not {"site", "color"}.issubset(reader.fieldnames or []):
            raise ValueError(f"{path} must contain columns: site, color")
        colors = {}
        for row in reader:
            site = (row.get("site") or "").strip()
            color = (row.get("color") or "").strip()
            if not site or not color:
                raise ValueError(f"{path}: site and color must be nonempty")
            if site in colors:
                raise ValueError(f"{path}: duplicate site {site}")
            colors[site] = color
        return colors


def parse_node_presence(path):
    if not path or not os.path.exists(path):
        return {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"node", "site", "observed", "mu"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"{path} must contain columns: {', '.join(sorted(required))}")
        presence = {}
        for row in reader:
            key = (clean_node_id(row["node"]), (row["site"] or "").strip())
            observed, selected = float(row["observed"]), float(row["mu"])
            if observed not in (0, 1) or selected not in (0, 1):
                raise ValueError(f"{path}: observed and mu must be binary")
            if key in presence:
                raise ValueError(f"{path}: duplicate node-site pair {key}")
            presence[key] = {"observed": int(observed), "mu": int(selected)}
        return presence


def parse_node_to_cluster(tree_path):
    if not tree_path or not os.path.exists(tree_path):
        return {}
    with open(tree_path) as fh:
        sections = fh.read().split("\n\n")
    node_to_cluster = {}
    for section in sections:
        lines = [line.strip() for line in section.splitlines() if line.strip()]
        if not lines:
            continue
        header = lines[0].split()
        upper = [h.upper() for h in header]
        if "CLUSTER" not in upper or "NODE_OF_CLUSTER" not in upper:
            continue
        cluster_idx = upper.index("CLUSTER")
        node_idx = upper.index("NODE_OF_CLUSTER")
        need = max(cluster_idx, node_idx)
        for line in lines[1:]:
            parts = line.split()
            if len(parts) <= need:
                continue
            node_to_cluster[clean_node_id(parts[node_idx])] = parts[cluster_idx]
    return node_to_cluster


def infer_tree_path(edge_path):
    base = os.path.basename(edge_path)
    match = EDGE_RE.match(base)
    if not match:
        return None
    return os.path.join(
        os.path.dirname(edge_path),
        f"{match.group('prefix')}.tree_{match.group('idx')}{match.group('suffix')}.txt",
    )


def infer_nodes_path(edge_path):
    match = EDGE_RE.match(os.path.basename(edge_path))
    if not match:
        return None
    return os.path.join(
        os.path.dirname(edge_path),
        f"{match.group('prefix')}.refined_nodes_sol{match.group('idx')}{match.group('suffix')}.tsv",
    )


def make_site_color_map(site_order, *site_sets, palette=None, overrides=None):
    ordered = []
    seen = set()
    for site in site_order or []:
        if site and site not in seen:
            seen.add(site)
            ordered.append(site)
    extras = set()
    for sites in site_sets:
        extras.update(s for s in sites if s)
    for site in sorted(extras - seen, key=natural_key):
        seen.add(site)
        ordered.append(site)
    palette = SITE_COLORS if palette is None else palette
    overrides = overrides or {}
    return {site: overrides.get(site, palette[i % len(palette)]) for i, site in enumerate(ordered)}


def migration_counts(edges):
    counts = Counter()
    for edge in edges:
        if edge.get("type") != "migration":
            continue
        source = edge.get("from_site")
        target = edge.get("to_site")
        if not source or not target or source == target:
            continue
        counts[source, target] += 1
    return counts


def write_counts_tsv(path, counts):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["from_site", "to_site", "count"])
        for (source, target), count in sorted(counts.items(), key=lambda x: natural_key(x[0][0]) + natural_key(x[0][1])):
            writer.writerow([source, target, count])


def render_dot(dot_path, pdf_path):
    subprocess.run(["dot", "-Tpdf", dot_path, "-o", pdf_path], check=True)


def write_dot(path, lines):
    with open(path, "w") as fh:
        fh.write("\n".join(lines))
        fh.write("\n")


def migration_node_key(node_id, site):
    return f"{node_id}|{site or 'NA'}"


def split_migration_node_key(key):
    if "|" not in str(key):
        return str(key), ""
    node_id, site = str(key).split("|", 1)
    return node_id, "" if site == "NA" else site


def migration_node_sort_key(key, node_to_cluster):
    node_id, site = split_migration_node_key(key)
    if str(node_id) == "0":
        cluster_part = [(-1, "")]
    else:
        cluster = node_to_cluster.get(node_id)
        cluster_part = [(0, part) for part in typed_natural_key(cluster if cluster is not None else node_id)]
    site_part = [(1, part) for part in typed_natural_key(site)]
    node_part = [(2, part) for part in typed_natural_key(node_id)]
    return cluster_part + site_part + node_part


def migration_edge_sort_key(edge, node_to_cluster):
    source = migration_node_key(edge["from_node"], edge["from_site"])
    target = migration_node_key(edge["to_node"], edge["to_site"])
    return migration_node_sort_key(source, node_to_cluster) + migration_node_sort_key(target, node_to_cluster)


def build_migration_dot(edges, node_to_cluster, site_color_map, title, node_presence=None,
                        style="clone-sites"):
    compact = style == "clone-sites"
    node_presence = node_presence or {}
    note = "Dotted node outline: site call changed from observed 0 to mu 1." if compact else None
    lines = [
        "digraph G {",
        "  graph [rankdir=TB, nodesep=0.35, ranksep=0.45, fontsize=18, labelloc=t, labeljust=l, bgcolor=white, ordering=out, label=" + graph_label_html(title, site_color_map) + "];",
        ('  node [shape=ellipse, style=filled, fillcolor=white, fontname=Helvetica, fontsize=14, '
         'width=0.55, height=0.35, margin="0.06,0.03", penwidth=2];') if compact else
        '  node [shape=box, style="rounded,filled", fontname=Helvetica, fontsize=14];',
        "  edge [fontname=Helvetica, arrowsize=0.7];" if compact else
        "  edge [fontname=Helvetica, arrowsize=0.8];",
    ]

    nodes = {}
    for edge in edges:
        for node_id, site in [(edge["from_node"], edge["from_site"]), (edge["to_node"], edge["to_site"])]:
            key = migration_node_key(node_id, site)
            if key in nodes:
                continue
            cluster = node_to_cluster.get(node_id)
            if compact:
                if cluster is not None:
                    label = str(cluster)
                elif str(node_id) == "0":
                    continue
                    label = "ROOT"
                else:
                    label = f"node {node_id}"
                border = site_color_map.get(site, "#555555")
                presence = node_presence.get((node_id, site), {})
                lifted = presence.get("observed") == 0 and presence.get("mu") == 1
                nodes[key] = {"label": label, "fillcolor": "white", "color": border,
                              "style": "filled,dashed" if lifted else "filled,solid",
                              "tooltip": f"Clone {cluster}, site {site}" if cluster is not None else label}
                continue
            if cluster is not None:
                label = f"C {cluster}\nS {site or 'NA'}"
            elif str(node_id) == "0":
                label = f"Root\nS {site or 'NA'}"
            else:
                label = f"C ?\nS {site or 'NA'}"
            fill = site_color_map.get(site, "lightyellow")
            nodes[key] = {"label": label, "fillcolor": fill}

    for key in sorted(nodes, key=lambda k: migration_node_sort_key(k, node_to_cluster)):
        lines.append(f"  {dot_quote(key)}{attr_text(nodes[key])};")

    for edge in sorted(edges, key=lambda e: migration_edge_sort_key(e, node_to_cluster)):
        if str(edge["from_node"]) == "0":
            continue
        source = migration_node_key(edge["from_node"], edge["from_site"])
        target = migration_node_key(edge["to_node"], edge["to_site"])
        etype = edge.get("type") or ""
        is_migration = etype == "migration"
        source_color = site_color_map.get(edge.get("from_site"), "black")
        target_color = site_color_map.get(edge.get("to_site"), source_color)
        if compact:
            source_color = nodes[source]["color"]
            target_color = nodes[target]["color"]
            color = (f"{source_color};0.5:{target_color};0.5"
                     if source_color != target_color else source_color)
            attrs = {"color": color, "penwidth": "1.8", "style": "solid"}
        elif is_migration:
            color = f"{source_color};0.5:{target_color};0.5"
            attrs = {"color": color, "penwidth": "2.5"}
        else:
            attrs = {"color": source_color, "penwidth": "1.6"}
        lines.append(f"  {dot_quote(source)} -> {dot_quote(target)}{attr_text(attrs)};")

    lines.append("}")
    return lines


def build_comigration_dot(counts, site_color_map, title):
    sites = set(site_color_map)
    for source, target in counts:
        sites.add(source)
        sites.add(target)
    max_count = max(counts.values(), default=1)
    lines = [
        "digraph G {",
        "  graph [rankdir=LR, nodesep=0.55, ranksep=0.75, fontsize=18, labelloc=t, labeljust=l, bgcolor=white, label=" + graph_label_html(title, site_color_map) + "];",
        "  node [shape=ellipse, style=filled, fontname=Helvetica, fontsize=14, penwidth=1.5];",
        "  edge [fontname=Helvetica, fontsize=12, arrowsize=0.9, color=firebrick, fontcolor=firebrick];",
    ]
    for site in sorted(sites, key=natural_key):
        attrs = {"label": site, "fillcolor": site_color_map.get(site, "lightyellow")}
        lines.append(f"  {dot_quote(site)}{attr_text(attrs)};")
    for (source, target), count in sorted(counts.items(), key=lambda x: natural_key(x[0][0]) + natural_key(x[0][1])):
        penwidth = 1.2 + (6.8 * float(count) / float(max_count))
        attrs = {"label": str(count), "penwidth": f"{penwidth:.2f}"}
        lines.append(f"  {dot_quote(source)} -> {dot_quote(target)}{attr_text(attrs)};")
    lines.append("}")
    return lines


def output_stem_for_edge(edge_path, output_dir, suffix):
    base = os.path.splitext(os.path.basename(edge_path))[0]
    return os.path.join(output_dir or os.path.dirname(edge_path), f"{base}.{suffix}")


def parse_edge_filename(edge_path):
    match = EDGE_RE.match(os.path.basename(edge_path))
    if not match:
        return None
    suffix = match.group("suffix")
    if suffix.endswith(".migration_graph") or suffix.endswith(".comigration_graph"):
        return None
    return match.group("prefix"), match.group("idx"), suffix


def discover_edge_files(path, include_best_migration=False):
    if os.path.isfile(path):
        files = [path]
    else:
        files = sorted(glob.glob(os.path.join(path, "*.refined_edges_sol*.tsv")), key=natural_key)
    if not include_best_migration:
        files = [f for f in files if "_best_migration" not in os.path.basename(f)]
    return [f for f in files if parse_edge_filename(f)]


def render_individual(edge_file, args, site_order):
    edges = parse_edges(edge_file)
    sites = {s for edge in edges for s in (edge.get("from_site"), edge.get("to_site")) if s}
    palette = CLONE_SITE_COLORS if args.migration_style == "clone-sites" else SITE_COLORS
    site_color_map = make_site_color_map(site_order, sites, palette=palette, overrides=args.site_colors)
    tree_path = infer_tree_path(edge_file)
    node_to_cluster = parse_node_to_cluster(tree_path)
    edge_base = os.path.splitext(os.path.basename(edge_file))[0]
    m = re.match(r'CRUK\d+', edge_base)
    if m:
        edge_base = m[0]

    if args.migration_graphs:
        nodes_path = infer_nodes_path(edge_file)
        node_presence = parse_node_presence(nodes_path) if args.migration_style == "clone-sites" else {}
        if args.migration_style == "clone-sites" and not node_presence:
            print(f"WARN: no CITUPm node presence sidecar for {edge_file}; cannot mark 0-to-1 site calls.", file=sys.stderr)
        stem = output_stem_for_edge(edge_file, args.output_dir, "migration_graph")
        dot_path = f"{stem}.dot"
        pdf_path = f"{stem}.pdf"
        write_dot(dot_path, build_migration_dot(edges, node_to_cluster, site_color_map, edge_base,
                                              node_presence=node_presence, style=args.migration_style))
        render_dot(dot_path, pdf_path)
        print(f"Saved migration graph to {pdf_path}")

    if args.comigration_graphs:
        counts = migration_counts(edges)
        stem = output_stem_for_edge(edge_file, args.output_dir, "comigration_graph")
        write_counts_tsv(f"{stem}.tsv", counts)
        dot_path = f"{stem}.dot"
        pdf_path = f"{stem}.pdf"
        write_dot(dot_path, build_comigration_dot(counts, site_color_map, edge_base))
        render_dot(dot_path, pdf_path)
        print(f"Saved comigration graph to {pdf_path}")

    return edges


def render_conjoined(edge_files, args, site_order):
    grouped = defaultdict(list)
    for path in edge_files:
        parsed = parse_edge_filename(path)
        if not parsed:
            continue
        prefix, idx, suffix = parsed
        key = (os.path.dirname(path), prefix, idx)
        grouped[key].append((suffix, path))

    for key in sorted(grouped, key=lambda k: natural_key(k[1]) + natural_key(k[2])):
        entries = sorted(grouped[key], key=lambda item: natural_key(item[0]))
        base = [(suffix, path) for suffix, path in entries if suffix == ""]
        refined = [(suffix, path) for suffix, path in entries if "_refined_" in suffix]
        selected = base + refined if refined else base
        if not selected:
            selected = entries

        total = Counter()
        sites = set(site_order)
        for _, path in selected:
            edges = parse_edges(path)
            total.update(migration_counts(edges))
            for edge in edges:
                for site in (edge.get("from_site"), edge.get("to_site")):
                    if site:
                        sites.add(site)

        directory, prefix, idx = key
        palette = CLONE_SITE_COLORS if args.migration_style == "clone-sites" else SITE_COLORS
        site_color_map = make_site_color_map(site_order, sites, palette=palette, overrides=args.site_colors)
        out_dir = args.output_dir or directory
        stem = os.path.join(out_dir, f"{prefix}.comigration_conjoined_sol{idx}")
        write_counts_tsv(f"{stem}.tsv", total)
        title = f"{prefix}.sol{idx} conjoined comigration ({len(selected)} refined graph" + ("" if len(selected) == 1 else "s") + ")"
        write_dot(f"{stem}.dot", build_comigration_dot(total, site_color_map, title))
        render_dot(f"{stem}.dot", f"{stem}.pdf")
        print(f"Saved conjoined comigration graph to {stem}.pdf")


def main():
    parser = argparse.ArgumentParser(description="Plot CITUPm migration and comigration graphs from CITUPm output.")
    parser.add_argument("path", help="CITUPm output directory or one *.refined_edges_sol*.tsv file")
    parser.add_argument("--sites-file", help="Optional sample/site TSV; site order controls colors")
    parser.add_argument("--site-colors", help="Optional TSV with site/color columns to override the palette")
    parser.add_argument("--migration-style", choices=["clone-sites", "detailed"], default="clone-sites",
                        help="Migration graph style: clone IDs in white site-colored ellipses (default), "
                             "or the previous detailed filled boxes.")
    parser.add_argument("--output-dir", help="Optional output directory")
    parser.add_argument("--include-best-migration-copy", action="store_true",
                        help="Also plot *_best_migration sidecar copies")
    parser.add_argument("--no-migration-graphs", dest="migration_graphs", action="store_false",
                        help="Skip expanded node-site migration graphs")
    parser.add_argument("--no-comigration-graphs", dest="comigration_graphs", action="store_false",
                        help="Skip per-file site-count comigration graphs")
    parser.add_argument("--no-conjoined", dest="conjoined", action="store_false",
                        help="Skip per-topology conjoined comigration graphs")
    parser.set_defaults(migration_graphs=True, comigration_graphs=True, conjoined=True)
    args = parser.parse_args()
    args.site_colors = parse_site_colors(args.site_colors)

    if not os.path.exists(args.path):
        parser.error(f"Path does not exist: {args.path}")
    if args.output_dir:
        os.makedirs(args.output_dir, exist_ok=True)

    edge_files = discover_edge_files(args.path, include_best_migration=args.include_best_migration_copy)
    if not edge_files:
        parser.error(f"No *.refined_edges_sol*.tsv files found under {args.path}")

    site_order = parse_site_order(args.sites_file)
    for edge_file in edge_files:
        render_individual(edge_file, args, site_order)
    if args.conjoined and os.path.isdir(args.path):
        render_conjoined(edge_files, args, site_order)


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as exc:
        sys.stderr.write(f"Graphviz dot failed with exit code {exc.returncode}: {' '.join(exc.cmd)}\n")
        sys.exit(exc.returncode)
