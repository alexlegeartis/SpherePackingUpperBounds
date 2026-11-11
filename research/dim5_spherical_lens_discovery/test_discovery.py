"""Portable deterministic smoke and rejection tests; no optimizer is run."""
if not __debug__:
    raise RuntimeError('Tests must not run under python -O.')

import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent
SUITE = ROOT.parent / 'spherical_lens_verification_package' if (ROOT.parent / 'spherical_lens_verification_package').is_dir() else ROOT.parents[1] / 'verifications/dim5/refinement_spherical_lens'
SOURCE_SHA256 = '28585450c46b4fbdc6f9eefec55b7e6d29f18fc527b654d00359ed425cb1f949'
CLI_FILES = ['reproduce_exact.py', 'run_search.py', 'export_and_check.py',
             'prepare_candidate.py', 'search_wedge_step.py', 'export_low.py']


class DiscoveryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='lens-discovery-test-')
        self.path = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def run_cli(self, name, *args, root=ROOT, optimized=False, success=True):
        result = subprocess.run([sys.executable, *(['-O'] if optimized else []),
                                 str(root / name), *map(str, args)],
                                cwd=self.path, text=True, capture_output=True)
        if success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def test_compile_and_cli_help(self):
        for path in ROOT.glob('*.py'):
            compile(path.read_bytes(), str(path), 'exec')
        for name in CLI_FILES:
            self.run_cli(name, '--help')

    def test_python_O_rejected(self):
        for path in ROOT.glob('*.py'):
            result = self.run_cli(path.name, '--help', optimized=True, success=False)
            self.assertIn('python -O', result.stderr)

    def test_exact_reproduction_and_no_overwrite(self):
        output = self.path / 'reproduced'
        self.run_cli('reproduce_exact.py', output)
        payload = (output / 'source_exact.json').read_bytes()
        self.assertEqual(hashlib.sha256(payload).hexdigest(), SOURCE_SHA256)
        self.assertEqual(payload, (ROOT / 'source_exact.json').read_bytes())
        self.assertTrue(all(r['off_origin_diagonal_shift'] == '0' for r in
                            json.loads(payload)[0]['exact_source_preparation']['repairs']))
        self.run_cli('reproduce_exact.py', output, success=False)
        self.assertEqual(payload, (output / 'source_exact.json').read_bytes())

    def test_final_search_parameters_and_local_paths(self):
        output = self.path / 'candidate.json'
        result = self.run_cli('run_search.py', '--dry-run', '--output', output)
        plan = json.loads(result.stdout)
        self.assertIs(plan['optimizer_run'], False)
        args = plan['command'][2:]
        for key, value in [('dimension', '5'), ('nodes', '5'), ('angular', '12'),
                           ('degree', '31'), ('ratio', '2.25'), ('subdivisions', '2'),
                           ('floor', '1e-9'), ('ridge', '0.03'), ('margin', '2'),
                           ('prox-floor', '0.0001'), ('rounds', '2')]:
            self.assertEqual(args[args.index('--' + key) + 1], value)
        for flag in ['--counts', '--schur', '--boundary']:
            self.assertIn(flag, args)
        for key, name in [('step-budget', 'step_allowances_spherical_005.json'),
                          ('far-boxes', 'far_boxes5_v2.json'), ('seed', 'wedge5_step.json')]:
            self.assertEqual(Path(args[args.index('--' + key) + 1]), ROOT / name)
        self.assertFalse(output.exists())
        self.run_cli('run_search.py', '--dry-run', '--rounds', '0', '--output', output, success=False)
        output.write_text('existing output\n')
        self.run_cli('run_search.py', '--dry-run', '--output', output, success=False)
        self.assertEqual(output.read_text(), 'existing output\n')

    def test_off_origin_shift_is_forbidden(self):
        source = json.loads((ROOT / 'wedge5_spherical.json').read_text())
        source[0]['matrices'][1][0][0] = -1
        candidate, output = self.path / 'invalid.json', self.path / 'prepared.json'
        candidate.write_text(json.dumps(source))
        result = self.run_cli('prepare_candidate.py', candidate, output, success=False)
        self.assertIn('off-origin shift', result.stderr)
        self.assertFalse(output.exists())

    def test_portable_canonical_export_without_optimization(self):
        relocated = self.path / 'portable'
        discovery = relocated / 'spherical_lens_discovery_package'
        proof = relocated / 'spherical_lens_verification_package'
        shutil.copytree(ROOT, discovery, ignore=shutil.ignore_patterns('__pycache__'))
        shutil.copytree(SUITE, proof, ignore=shutil.ignore_patterns('__pycache__'))
        output = self.path / 'exported'
        self.run_cli('export_and_check.py', '--export-only', '--output', output, root=discovery)
        report = json.loads((output / 'export_report.json').read_text())
        self.assertTrue(all(report['canonical_export_equality'].values()))
        self.assertIs(report['complete_verification_passed'], False)
        self.assertEqual((output / 'source_exact.json').read_bytes(), (ROOT / 'source_exact.json').read_bytes())
        self.run_cli('export_and_check.py', '--export-only', '--output', output,
                     root=discovery, success=False)
        # Reject mutated budget bytes before any output directory is created.
        (discovery / 'step_allowances_spherical_005.json').write_bytes(b'{}\n')
        result = self.run_cli('export_and_check.py', '--export-only',
                              '--output', self.path / 'bad-budget', root=discovery, success=False)
        self.assertIn('Wrong spherical lens budget', result.stderr)
        self.assertFalse((self.path / 'bad-budget').exists())
        result = self.run_cli('run_search.py', '--dry-run', '--output',
                              self.path / 'bad-search.json', root=discovery, success=False)
        self.assertIn('Wrong spherical lens budget', result.stderr)
        self.assertFalse((self.path / 'bad-search.json').exists())

    def test_fresh_candidate_is_not_automatically_certified(self):
        output = self.path / 'not-a-proof'
        result = self.run_cli('export_and_check.py', '--candidate', ROOT / 'wedge5_spherical.json',
                              '--output', output, success=False)
        self.assertIn('fresh candidates are not certified', result.stderr)
        self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main(verbosity=2)
