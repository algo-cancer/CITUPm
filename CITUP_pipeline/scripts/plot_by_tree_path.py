import argparse
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Plot each person's trees in a cohort directory.",
    )
    parser.add_argument(
        "directory",
        type=Path,
        help="Cohort directory containing PERSON/trees/ subdirectories.",
    )
    args = parser.parse_args()

    directory = args.directory.resolve()
    if not directory.is_dir():
        parser.error(f"Not a directory: {args.directory}")

    cohort = directory.name

    for child in sorted(directory.iterdir()):
        if not child.is_dir():
            continue

        indir = child / "trees"
        if not indir.is_dir():
            print(f"Skipping {child}: no trees directory", file=sys.stderr)
            continue

        person = child.name
        outdir = child / "plots"
        sites = Path("data") / cohort / person / "sites.tsv"

        subprocess.run(
            [
                sys.executable,
                "../plot_with_gene_names.py",
                str(indir),
                "--output-dir",
                str(outdir),
                "--sample-sites",
                str(sites),
                "--with-site-transitions",
                "--with_node_info",
                "--all-migration-dag-colorings",
            ],
            check=True,
        )
        if any(indir.glob("*.refined_edges_sol*.tsv")):
            subprocess.run(
                [
                    sys.executable,
                    "../plot_comigration_graphs.py",
                    str(indir),
                    "--output-dir",
                    str(outdir),
                    "--sites-file",
                    str(sites),
                ],
                check=True,
            )


if __name__ == "__main__":
    main()
