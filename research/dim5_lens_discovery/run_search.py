"""Run only the bounded final 5D warm-start search; numerical output is not a bound."""
if not __debug__:
    raise RuntimeError('Discovery safety checks must not run under python -O.')

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', type=Path, default=ROOT / 'wedge5.json')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--rounds', type=int, default=2)
    parser.add_argument('--dry-run', action='store_true',
                        help='Print the complete command without running an optimizer')
    args = parser.parse_args()
    output, seed = args.output.resolve(), args.seed.resolve()
    if output.exists() or args.rounds < 1:
        parser.error('Use a new output path and a positive round count')
    if not seed.is_file():
        parser.error('Seed file does not exist')
    for name in ['cyclic_fiber_rows_5.npz', 'far_boxes5_v2.json', 'step_allowances_005.json']:
        if not (ROOT / name).is_file():
            parser.error('Missing bundled input: ' + name)
    budget_hash = hashlib.sha256((ROOT / 'step_allowances_005.json').read_bytes()).hexdigest()
    if budget_hash != 'cd4af8aaac2b6d2b1f0559a843953ee0a9132543fb371a3fad4ab9d14e0a44e2':
        parser.error('Wrong centroid lens budget; this final-stage search uses the frozen old caps')
    command = [sys.executable, str(ROOT / 'search_wedge_step.py'),
               '--dimension', '5', '--nodes', '5', '--angular', '12',
               '--degree', '31', '--ratio', '2.25', '--subdivisions', '2',
               '--counts', '--schur', '--floor', '1e-9', '--ridge', '0.2',
               '--margin', '2', '--prox-floor', '0.0001', '--boundary',
               '--exclusion', '1.0', '--wedge-radius', '1.08',
               '--far-boxes', str(ROOT / 'far_boxes5_v2.json'),
               '--step-budget', str(ROOT / 'step_allowances_005.json'),
               '--rounds', str(args.rounds), '--seed', str(seed), '--output', str(output)]
    if args.dry_run:
        print(json.dumps({'command': command, 'optimizer_run': False,
                          'warning': 'Numerical candidates are not certified bounds.'}, indent=2))
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(command, check=True)
    print('NUMERICAL_ONLY: exact preparation and complete verification remain mandatory.')


if __name__ == '__main__':
    main()
