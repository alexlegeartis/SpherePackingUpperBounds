"""Re-run one successful numerical stage. Output is never a certificate."""
import argparse
from pathlib import Path
from search_joint import solve

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('stage', choices=['joint', 'refine'])
ap.add_argument('--seed', type=Path)
ap.add_argument('--output', type=Path, required=True)
ap.add_argument('--rounds', type=int)
args = ap.parse_args()
root = Path(__file__).resolve().parent
joint = args.stage == 'joint'
seed = args.seed or root / ('candidate_certificate_v2_6.json' if joint else 'seed_feasible5.json')
rounds = args.rounds if args.rounds is not None else (3 if joint else 4)
if rounds < 1 or args.output.exists():
    ap.error('Use a positive round count and a new output path')
args.output.parent.mkdir(parents=True, exist_ok=True)
solve(n=6, degree=31, nodes=5, angular=12, ratio=2.25, subdivisions=2,
      rounds=rounds, pairwise=False, output=str(args.output.resolve()),
      exclusion=1.0, counts=True, seed=str(seed.resolve()), floor=1e-8,
      schur=True, joint=joint, ridge=0.0001 if joint else 0.01,
      margin=1.0 if joint else 3.0)
print('NUMERICAL_ONLY: complete exact verification is mandatory.')
