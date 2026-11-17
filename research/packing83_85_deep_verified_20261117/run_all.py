from pathlib import Path
import argparse,hashlib,json,shutil,tempfile,subprocess,sys
p=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--full',action='store_true');args=ap.parse_args()
if not __debug__:raise RuntimeError('Assertions required')
for fn,h in json.loads((p/'SHA256SUMS.json').read_text()).items():assert hashlib.sha256((p/fn).read_bytes()).hexdigest()==h,fn
with tempfile.TemporaryDirectory(prefix='packing82plane_')as tmp:
 w=Path(tmp)/'work';shutil.copytree(p,w)
 def run(cmd):print('RUN',*map(str,cmd),flush=True);subprocess.run(list(map(str,cmd)),check=True,cwd=w)
 run([sys.executable,w/'run82.py']+(['--full']if args.full else[]))
 for f in [w/'l80_diagonal_automorphisms_20261117/verify_plane.py',w/'l80_diagonal_automorphisms_20261117/build82_plane.py',w/'l80_plane82_independent_20261117/check.py']:run([sys.executable,f])
 print('PASS improved82: independent full symmetry transfers, XOR plane, A2 Gram/parity, exact density floor, inherited complete coset proofs')
