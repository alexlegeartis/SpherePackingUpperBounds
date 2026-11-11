"""Export prerequisites, optionally replaying a supplied frozen verifier suite.

All output goes to a new directory. No source verifier or source candidate is
modified. Without --verification-package this is NOT a complete verification.
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parent
    parser.add_argument('--candidate', type=Path, default=root / 'mix_continue_99hundredths.json')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--verification-package', type=Path)
    args = parser.parse_args()
    candidate = args.candidate.resolve()
    output = args.output.resolve()
    if output.exists():
        parser.error('Output already exists; choose a new directory')
    if not candidate.is_file():
        parser.error('Candidate file does not exist')
    verifier_files = ['verify.py', 'gaussian_guards.py',
                      'certify_count_pair_arb.py', 'certify_far_bernstein.py',
                      'certify_cyclic_bernstein.py']
    suite = args.verification_package.resolve() if args.verification_package else None
    if suite and any(not (suite / f).is_file() for f in verifier_files):
        parser.error('The supplied verification package is incomplete')
    digest = hashlib.sha256(candidate.read_bytes()).hexdigest()
    output.mkdir(parents=True)
    subprocess.run([sys.executable, str(root / 'export_cheb.py'),
                    str(candidate), str(output)], check=True)
    if hashlib.sha256(candidate.read_bytes()).hexdigest() != digest:
        raise RuntimeError('Candidate changed during export')
    report = {'candidate_sha256': digest, 'prerequisite_export_passed': True,
              'complete_verification_passed': False}
    if suite:
        for name in verifier_files:
            source_name = 'verify_primary.py' if name == 'verify.py' and (suite / 'verify_primary.py').is_file() else name
            shutil.copyfile(suite / source_name, output / name)
        # The frozen 6D driver also enforces its strict upper decimal.
        subprocess.run([sys.executable, str(output / 'verify.py')], check=True)
        result = json.loads((output / 'verification_result.json').read_text())
        if result.get('verified') is not True:
            raise RuntimeError('Complete verification did not pass')
        report['complete_verification_passed'] = True
        report['result'] = result
    (output / 'export_report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    if not suite:
        print('PREREQUISITES_ONLY: complete pair/far/cyclic verification is still required.')


if __name__ == '__main__':
    main()
