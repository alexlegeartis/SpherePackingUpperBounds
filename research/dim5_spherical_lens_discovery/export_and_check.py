"""Re-export and verify the frozen 5D source using its sibling proof package.

The verifier is deliberately bound to the frozen source hash. This command
does not certify arbitrary fresh numerical outputs. All generated files go
to a new directory; the proof package and discovery inputs are never edited.
"""
if not __debug__:
    raise RuntimeError('Verification must not run under python -O.')

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
SOURCE_SHA256 = '28585450c46b4fbdc6f9eefec55b7e6d29f18fc527b654d00359ed425cb1f949'
BUDGET_SHA256 = '605a15c07f305a590bb1e0e2c80a9a88f3991023129bcbb5a88c103aa4552046'
VERIFIER_FILES = ['verify.py', 'gaussian_guards.py', 'audit_algebra_exact.py',
                  'construct_step_allowances.py', 'step_allowance_query.py',
                  'verify_count_endpoints.py', 'certify_count_pair_arb.py',
                  'certify_cyclic_bernstein.py', 'certify_general_lens_far.py', 'spherical_lens_caps.py']
EXPORT_FILES = ['witness.json', 'witness_psd_exact.json', 'prerequisites.json']


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate', type=Path, default=ROOT / 'source_exact.json',
                        help='Frozen exact source or its byte-identical reconstruction')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--verification-package', type=Path,
                        default=(ROOT.parent / 'spherical_lens_verification_package' if (ROOT.parent / 'spherical_lens_verification_package').is_dir() else ROOT.parents[1] / 'verifications/dim5/refinement_spherical_lens'))
    parser.add_argument('--export-only', action='store_true',
                        help='Check canonical export equality without complete verification')
    args = parser.parse_args()
    candidate, output, suite = (p.resolve() for p in
                                (args.candidate, args.output, args.verification_package))
    if output.exists():
        parser.error('Output already exists; choose a new directory')
    required = VERIFIER_FILES + EXPORT_FILES + ['source_exact.json', 'step_allowances_spherical_005.json']
    if any(not (suite / name).is_file() for name in required):
        parser.error('The supplied lens verification package is incomplete')
    if not candidate.is_file() or digest(candidate) != SOURCE_SHA256:
        parser.error('This replay requires the frozen exact source; fresh candidates are not certified')
    if digest(suite / 'source_exact.json') != SOURCE_SHA256:
        raise RuntimeError('Sibling verifier source hash mismatch')
    for budget in [ROOT / 'step_allowances_spherical_005.json', suite / 'step_allowances_spherical_005.json']:
        if digest(budget) != BUDGET_SHA256:
            raise RuntimeError('Wrong spherical lens budget: ' + str(budget))
    output.mkdir(parents=True)
    subprocess.run([sys.executable, str(ROOT / 'export_low.py'),
                    str(candidate), str(output)], check=True)
    # Preserve all mathematical fields. Only the lens metadata differs from
    # the low-dimensional baseline exporter, matching the approved producer.
    witness = json.loads((output / 'witness.json').read_text())
    witness.update(kernel_family='tail_funded_lens_cyclic',
                   scope='Research candidate using geometric lens counts; requires all checks and independent review.',
                   lens_budget='step_allowances_spherical_005.json')
    (output / 'witness.json').write_text(json.dumps(witness, indent=2) + '\n')
    canonical_matches = {}
    for name in EXPORT_FILES:
        canonical_matches[name] = (json.loads((output / name).read_text()) ==
                                   json.loads((suite / name).read_text()))
        if not canonical_matches[name]:
            raise RuntimeError('Canonical export mismatch: ' + name)
    shutil.copyfile(candidate, output / 'source_exact.json')
    shutil.copyfile(ROOT / 'step_allowances_spherical_005.json', output / 'step_allowances_spherical_005.json')
    for name in VERIFIER_FILES:
        source_name = 'verify_primary.py' if name == 'verify.py' and (suite / 'verify_primary.py').is_file() else name
        shutil.copyfile(suite / source_name, output / name)
    if digest(candidate) != SOURCE_SHA256:
        raise RuntimeError('Candidate changed during export')
    report = {'source_sha256': SOURCE_SHA256, 'spherical_budget_sha256': BUDGET_SHA256,
              'canonical_export_equality': canonical_matches,
              'complete_verification_passed': False,
              'verifier_script_sha256': {name: digest(output / name) for name in VERIFIER_FILES}}
    if not args.export_only:
        subprocess.run([sys.executable, str(output / 'verify.py')], check=True, cwd=output)
        result = json.loads((output / 'verification_result.json').read_text())
        if result.get('verified') is not True or result.get('strict_upper') != '0.509953498':
            raise RuntimeError('Complete verification did not pass the frozen strict upper bound')
        report.update(complete_verification_passed=True, result=result)
    (output / 'export_report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    if args.export_only:
        print('PREREQUISITES_ONLY: complete verification was not run.')
    else:
        print('FROZEN_EXPORT_AND_COMPLETE_VERIFICATION_PASS')


if __name__ == '__main__':
    main()
