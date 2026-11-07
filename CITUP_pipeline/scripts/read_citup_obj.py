import sys
import re
import hashlib

cruk = re.compile(r"CRUKP?\d+(?:_Tumour\d)?")

for filename in sys.argv[1:]:
    m = cruk.search(filename)
    if not m:
        print("Cannot read", filename)
        sys.exit(0)
    person = m.group()

    score = "NA"
    migrations = "NA"
    comigrations = "NA"
    tree = []
    reading_score = False
    reading_migrations = False
    reading_comigrations = False
    reading_tree = False
    with open(filename) as fh:
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
                if stripped_line == '':
                    reading_tree = False
                    continue
                tree.append(stripped_line)
            elif line.startswith("TREE_SCORE_NO_REWARD"):
                reading_score = True
            elif line.startswith("MIKE_MIGRATION_SCORE"):
                reading_migrations = True
            elif line.startswith("MIKE_COMIGRATION_SCORE"):
                reading_comigrations = True
            elif line.startswith('PARENT_NODE'):
                reading_tree = True
            else:
                reading_score = False
                reading_migrations = False
                reading_tree = False
    tree_digest = hashlib.md5('\n'.join(tree).encode("utf-8")).hexdigest()

    print(person, migrations, comigrations, score, tree_digest, sep="\t")
