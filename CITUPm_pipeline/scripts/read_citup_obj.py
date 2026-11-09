import argparse
import hashlib
from pathlib import Path


parser = argparse.ArgumentParser(
    usage="%(prog)s DIRECTORY",  # Edit this usage line as needed.
    description="Extract scores and tree digests from *.tree_0.txt files.",
)
parser.add_argument(
    "directory",
    type=Path,
    help="Directory to search recursively.",
)
args = parser.parse_args()


for filename in args.directory.rglob("*.tree_0.txt"):
    person = filename.name[:-len(".tree_0.txt")]

    score = "NA"
    migrations = "NA"
    comigrations = "NA"
    tree = []
    reading_score = False
    reading_migrations = False
    reading_comigrations = False
    reading_tree = False

    with filename.open() as fh:
        for line in fh:
            if reading_score:
                score = line.rstrip()
                reading_score = False
            elif reading_migrations:
                migrations = line.rstrip()
                reading_migrations = False
            elif reading_comigrations:
                comigrations = line.rstrip()
                reading_comigrations = False
            elif reading_tree:
                stripped_line = line.rstrip()
                if stripped_line == "":
                    reading_tree = False
                    continue
                tree.append(stripped_line)
            elif line.startswith("TREE_SCORE_NO_REWARD"):
                reading_score = True
            elif line.startswith("NUM_MIGRATIONS"):
                reading_migrations = True
            elif line.startswith("NUM_COMIGRATIONS"):
                reading_comigrations = True
            elif line.startswith("PARENT_NODE"):
                reading_tree = True
            else:
                reading_score = False
                reading_migrations = False
                reading_tree = False

    tree_digest = hashlib.md5(
        "\n".join(tree).encode("utf-8")
    ).hexdigest()

    print(person, migrations, comigrations, score, tree_digest, sep="\t")
