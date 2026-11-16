"""Materialize hash-checked certificates in a temporary copy and replay them."""
import argparse,gzip,hashlib,json,shutil,subprocess,sys,tempfile
from pathlib import Path
if not __debug__:raise RuntimeError('Assertions required; do not run with -O')
p=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--dimension',type=int,choices=[14,15,16,17,18,19,20,21,22,23,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48]);a=ap.parse_args()
m=json.loads((p/'MANIFEST.json').read_text())
with tempfile.TemporaryDirectory(prefix='lattice-proof-') as tmp:
 root=Path(tmp)
 for entry in m['files']:
  rel=Path(entry['stored_path']);dest=Path(entry['path'])
  if rel.is_absolute() or dest.is_absolute() or '..' in rel.parts or '..' in dest.parts:raise ValueError('Unsafe manifest path')
  raw=(p/rel).read_bytes()
  if hashlib.sha256(raw).hexdigest()!=entry['stored_sha256']:raise ValueError('Stored hash mismatch: '+str(rel))
  data=gzip.decompress(raw) if entry['gzip'] else raw
  if hashlib.sha256(data).hexdigest()!=entry['sha256']:raise ValueError('Decompressed hash mismatch: '+str(dest))
  target=root/dest;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
 dims=[a.dimension] if a.dimension else [14,15,16,17,18,19,20,21,22,23,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48]
 for n in dims:
  print('REPLAY lattice dimension',n,flush=True)
  folder='dim48_cyclic' if n==48 else 'dim42_own' if n==42 else f'dim{n}_tight' if n in [38,39,40] else f'dim{n}_threewindow' if n in [32,35,36] else f'dim{n}'
  subprocess.run([sys.executable,str(root/folder/'run_all.py')],cwd=root/folder,check=True)
 print('ALL REQUESTED LATTICE REPLAYS PASSED',flush=True)
