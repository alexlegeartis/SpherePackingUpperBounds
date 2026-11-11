"""Run the direct follow-on discovery stage; output is not a certificate."""
import argparse
from pathlib import Path
import subprocess
import sys
from search_joint import solve

root = Path(__file__).resolve().parent
ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('--seed', type=Path, default=root / 'continue5_round3_snapshot.json')
ap.add_argument('--output', type=Path, required=True)
ap.add_argument('--rounds', type=int, default=3)
args = ap.parse_args()
if args.output.exists() or args.rounds < 1:
    ap.error('Use a new output path and positive round count')
if not (root / 'cyclic_fiber_rows.npz').exists():
    subprocess.run([sys.executable, str(root / 'build_cyclic_fibers.py')], check=True)
args.output.parent.mkdir(parents=True, exist_ok=True)
solve(n=6, degree=31, nodes=5, angular=12, ratio=2.25, subdivisions=2,
      rounds=args.rounds, pairwise=False, output=str(args.output.resolve()),
      exclusion=1.0, counts=True, seed=str(args.seed.resolve()), floor=1e-8,
      schur=True, joint=False, ridge=1.0, margin=3.0, boundary=True,
      far_patch=True, prox_floor=0.0001)
print('NUMERICAL_ONLY: complete exact verification is mandatory.')
