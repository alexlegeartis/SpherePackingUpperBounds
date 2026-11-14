"""Replay the sixteen registered profile packages, including registry/hash gates."""
import argparse,concurrent.futures,subprocess,sys
from pathlib import Path
if not __debug__:raise RuntimeError('Assertions must remain enabled')
ROOT=Path(__file__).resolve().parents[2]
ap=argparse.ArgumentParser();ap.add_argument('--jobs',type=int,default=1);args=ap.parse_args()
assert 1<=args.jobs<=4
dimensions=[22,25,32,*range(35,48)]
def replay(n):
    directory=ROOT/f'verifications/dim{n}/refinement_profile20261114'
    result=subprocess.run([sys.executable,str(directory/'verify.py')],cwd=directory,capture_output=True,text=True)
    print(result.stdout,flush=True)
    assert result.returncode==0,(n,result.stderr)
with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:list(pool.map(replay,dimensions))
print('PASS: all sixteen registered profile packages',flush=True)
