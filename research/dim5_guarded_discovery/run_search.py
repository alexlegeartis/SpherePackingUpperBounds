"""Run the bounded 5D final warm-start search; output is not a certificate."""
import argparse
from pathlib import Path
import subprocess,sys
from search_low_generic import solve
root=Path(__file__).resolve().parent
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--seed',type=Path,default=root/'seed_patch_dim5.json')
ap.add_argument('--output',type=Path,required=True)
ap.add_argument('--rounds',type=int,default=3)
args=ap.parse_args()
if args.output.exists() or args.rounds<1: ap.error('Use a new output path and positive round count')
if not (root/'cyclic_fiber_rows_5.npz').exists():
 subprocess.run([sys.executable,str(root/'build_low_fibers.py'),'--dimension','5'],check=True)
args.output.parent.mkdir(parents=True,exist_ok=True)
solve(n=5,degree=31,nodes=5,angular=12,ratio=2.25,subdivisions=2,
 rounds=args.rounds,pairwise=False,output=str(args.output.resolve()),exclusion=1.0,
 counts=True,seed=str(args.seed.resolve()),floor=1e-8,schur=True,joint=False,
 ridge=1.0,margin=3.0,prox_floor=0.0001,boundary=True,far_boxes=str(root/'far_boxes5.json'))
print('NUMERICAL_ONLY: complete exact verification remains mandatory.')
