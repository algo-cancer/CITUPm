"""Small end-to-end checks of the public CLI and emitted migration histories."""
import csv
import importlib.util
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


def read_header(path):
    lines = path.read_text().split('\n\n')[0].splitlines()
    return dict(zip(lines[::2], lines[1::2]))


def read_tsv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


@unittest.skipUnless(importlib.util.find_spec('gurobipy'), 'gurobipy is required')
class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='citupm-test-')
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)

    def run_cli(self, name, *options, frequencies=None):
        prefix = self.work / name / 'result'
        prefix.parent.mkdir(parents=True)
        cmd = [sys.executable, str(ROOT / 'citupm.py'),
               '-f', str(frequencies or ROOT / 'examples/frequencies.tsv'),
               '-o', str(prefix), '--objective', 'cluster', '--clusterWeight', 'equal',
               '--presence-threshold', '0.025', '--numSolutions', '1',
               '--gurobi-threads', '1', *options]
        run = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        self.assertEqual(run.returncode, 0, run.stdout[-3000:] + run.stderr)
        self.assertNotIn('<gurobi.', run.stdout)
        self.assertNotIn("[{'core':", run.stdout)
        return prefix

    def assert_history(self, prefix, migrations, comigrations, suffix=''):
        header = read_header(Path(f'{prefix}.tree_0{suffix}.txt'))
        self.assertEqual(int(header['NUM_MIGRATIONS']), migrations)
        self.assertEqual(int(header['NUM_COMIGRATIONS']), comigrations)
        self.assertEqual(header['MIGRATION_STATUS'], 'ok')
        self.assertEqual(header['SITE_DAG_FINAL'], 'YES')
        for key in ('MIGRATION_SCORE_SOURCE', 'MIGRATION_SCORE_REFERENCE',
                    'MIGRATION_SCORE_UPPER_BOUND', 'MIGRATION_BOUNDED_SCORE',
                    'SUM_M_MINUS_I', 'MU_ASSIGNMENT_INDEX', 'MU_ASSIGNMENT_COUNT',
                    'MU_ASSIGNMENT_DEDUP_REMOVED'):
            self.assertNotIn(key, header)
        self.assertFalse(any('MIKE' in key for key in header))
        edges = read_tsv(Path(f'{prefix}.refined_edges_sol0{suffix}.tsv'))
        calls = read_tsv(Path(f'{prefix}.refined_nodes_sol0{suffix}.tsv'))
        migration_edges = [e for e in edges if e['type'] == 'migration']
        self.assertEqual(len(migration_edges), migrations)
        self.assertEqual(len({(e['from_site'], e['to_site']) for e in migration_edges}), comigrations)
        self.assertTrue(all(e['from_node'] == e['to_node'] != '0' for e in migration_edges))
        self.assertTrue(all(e['from_site'] == '' for e in edges if e['from_node'] == '0'))
        self.assertEqual(sum(int(row['mu']) for row in calls) - 2, migrations)
        self.assertTrue(all(row['node'] != '0' for row in calls))
        return header, calls

    def test_paper_three_stage_workflow_and_plots(self):
        fit = self.run_cli('fit')
        fit_header = read_header(Path(f'{fit}.tree_0.txt'))
        self.assertEqual(fit_header['FIT_OPTIMAL'], 'YES')
        self.assertEqual(fit_header['NUM_MIGRATIONS'], 'NA')
        reference = float(fit_header['TREE_SCORE_NO_REWARD'])
        self.assertLess(abs(reference), 1e-8)
        joint = self.run_cli('joint', '--sites-file', str(ROOT / 'examples/sites.tsv'),
                             '--metastatic-sites', 'metastasis',
                             '--minimize-migrations-at-score', str(max(0, reference)),
                             '--max-refined-solutions', '0')
        header, _ = self.assert_history(joint, 2, 1)
        self.assertEqual(header['MIGRATION_MINIMUM_PROVEN'], 'YES')
        self.assertEqual(header['MIGRATION_SCORE_TIEBREAK_PROVEN'], 'YES')
        self.assertEqual(header['REFINED_ENUMERATION_EXHAUSTIVE'], 'YES')
        self.assertLessEqual(float(header['TREE_SCORE_NO_REWARD']),
                             (math.sqrt(max(0, reference)) + 0.05)**2 + 1e-8)
        scores = read_tsv(Path(f'{joint}.treeScores.tsv'))
        self.assertEqual(scores[0]['numMigrations'], '2')
        self.assertEqual(scores[0]['numComigrations'], '1')
        audit = dict((row['key'], row['value']) for row in
                     read_tsv(Path(f'{joint}.migration_score_optimization.tsv')))
        self.assertEqual(audit['migration_count'], '2')
        self.assertAlmostEqual(float(audit['score_upper_bound']),
                               (math.sqrt(max(0, reference)) + 0.05)**2)
        if shutil.which('dot') and all(importlib.util.find_spec(m) for m in ('graphviz', 'pandas')):
            for script, args in (
                ('plot_with_gene_names.py', ['--sample-sites', str(ROOT / 'examples/sites.tsv'),
                                             '--with-site-transitions', '--with-node-info']),
                ('plot_comigration_graphs.py', ['--sites-file', str(ROOT / 'examples/sites.tsv')]),
            ):
                run = subprocess.run([sys.executable, str(ROOT / script), str(joint.parent), *args],
                                     capture_output=True, text=True, timeout=60)
                self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            for path in (Path(f'{joint}.tree_0.pdf'),
                         Path(f'{joint}.refined_edges_sol0.migration_graph.pdf'),
                         Path(f'{joint}.refined_edges_sol0.comigration_graph.pdf'),
                         Path(f'{joint}.comigration_conjoined_sol0.pdf')):
                self.assertTrue(path.read_bytes().startswith(b'%PDF'), str(path))
            # The same tree plotter can consume a renamed edge sidecar directly.
            direct_pdf = self.work / 'direct.pdf'
            run = subprocess.run([sys.executable, str(ROOT / 'plot_with_gene_names.py'),
                                  f'{joint}.refined_edges_sol0.tsv', '--output', str(direct_pdf)],
                                 capture_output=True, text=True, timeout=60)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            self.assertTrue(direct_pdf.read_bytes().startswith(b'%PDF'))

    def test_recursive_admissibility_supports_unobserved_ancestor(self):
        frequencies = self.work / 'recursive.tsv'
        frequencies.write_text('mutID\tprimary_1\tprimary_2\tmetastasis_1\tclusterID\n'
                               'm1\t0.8\t0.8\t0.4\tC1\n'
                               'm2\t0.4\t0.4\t0.4\tC2\n')
        prefix = self.run_cli('recursive', '--sites-file', str(ROOT / 'examples/sites.tsv'),
                              '--metastatic-sites', 'metastasis', '--num-migrations', '2',
                              frequencies=frequencies)
        _, calls = self.assert_history(prefix, 2, 1)
        ancestor = next(row for row in calls if row['node'] == '1' and row['site'] == 'metastasis')
        self.assertEqual(ancestor['observed'], '0')
        self.assertEqual(ancestor['mu'], '1')
        self.assertAlmostEqual(float(ancestor['presence_signal']), 0)

    def test_all_minimum_comigration_refinements(self):
        prefix = self.run_cli('refinements', '--sites-file', str(ROOT / 'examples/sites.tsv'),
                              '--max-migrations', '2', '--max-refined-solutions', '0')
        header, _ = self.assert_history(prefix, 2, 1)
        self.assertEqual(header['REFINED_SOLUTION_COUNT'], '2')
        self.assertEqual(header['REFINED_ENUMERATION_EXHAUSTIVE'], 'YES')
        self.assert_history(prefix, 2, 1, '_refined_1')
        summaries = read_tsv(Path(f'{prefix}.refined_solution_summary_sol0.tsv'))
        self.assertEqual(len(summaries), 2)
        self.assertTrue(all(row['num_migrations'] == '2' and row['num_comigrations'] == '1'
                            for row in summaries))
        arcs = {row['location_graph_edges'] for row in summaries}
        self.assertEqual(arcs, {'primary->metastasis', 'metastasis->primary'})

    def test_skipping_stage_three_does_not_claim_fit_optimality(self):
        prefix = self.run_cli('migration_only', '--sites-file', str(ROOT / 'examples/sites.tsv'),
                              '--metastatic-sites', 'metastasis',
                              '--minimize-migrations-at-score', '0',
                              '--no-migration-score-tiebreak')
        header, _ = self.assert_history(prefix, 2, 1)
        self.assertEqual(header['MIGRATION_MINIMUM_PROVEN'], 'YES')
        self.assertEqual(header['MIGRATION_SCORE_TIEBREAK_ENABLED'], 'NO')
        self.assertEqual(header['MIGRATION_SCORE_TIEBREAK_PROVEN'], 'NO')
        self.assertEqual(header['FIT_OPTIMAL'], 'NO')

    def test_distinct_topology_enumeration(self):
        frequencies = self.work / 'ties.tsv'
        frequencies.write_text('mutID\tprimary_1\tprimary_2\tmetastasis_1\tclusterID\n'
                               'm1\t0.6\t0.6\t0.5\tC1\n'
                               'm2\t0.6\t0.6\t0.5\tC2\n')
        prefix = self.run_cli('ties', '--numSolutions', '2', frequencies=frequencies)
        scores = read_tsv(Path(f'{prefix}.treeScores.tsv'))
        self.assertEqual(len(scores), 2)
        sections = [Path(f'{prefix}.tree_{i}.txt').read_text().split('\n\n')[1]
                    for i in range(2)]
        self.assertNotEqual(*sections)
        self.assertTrue(all(abs(float(row['fitWithoutReward'])) < 1e-6 for row in scores))


if __name__ == '__main__':
    unittest.main()
