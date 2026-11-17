"""Portable exact/interval verification. Assertions must be enabled."""
import argparse,json,hashlib,sys,subprocess,tempfile,shutil,gzip
from pathlib import Path
if not __debug__:raise RuntimeError('Do not use python -O')
D=Path(__file__).resolve().parent;ap=argparse.ArgumentParser();ap.add_argument('--quick',action='store_true',help='Exact algebra/cap only; NOT a full proof replay');ap.add_argument('--dimension',type=int,choices=[45,46,47,48]);a=ap.parse_args();dims=[a.dimension]if a.dimension else[45,46,47,48]
manifest=json.loads((D/'MANIFEST.json').read_text());stored=json.loads((D/'STORED_FILES.json').read_text())
with tempfile.TemporaryDirectory(prefix='compact-power6-proof-')as td:
 P=Path(td)/'proof';P.mkdir()
 for f,h in manifest.items():
  entry=stored.get(f,{});sf=entry.get('stored_path',f)
  for rel in [f,sf]:
   rp=Path(rel)
   if rp.is_absolute() or '..' in rp.parts:raise ValueError('Unsafe manifest path')
  data=(D/sf).read_bytes()
  if entry:
   assert hashlib.sha256(data).hexdigest()==entry['stored_sha256'],('stored hash mismatch',f)
   data=gzip.decompress(data)
  assert hashlib.sha256(data).hexdigest()==h,('raw hash mismatch',f)
  target=P/f;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
 def run(f,*args):
  r=subprocess.run([sys.executable,str(P/f),*map(str,args)],cwd=P,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);print(r.stdout,flush=True);assert r.returncode==0,('failed',f)
 run('source48/check_exact.py');run('source48/cap_independent.py');run('transfers/check_exact.py')
 if not a.quick:
  assert (P/'source48/input/main_result.json').exists(),'Main geometry incomplete; no full certificate yet'
  run('source48/replay_trace.py','main')
  for n in sorted(set(dims+[48])):
   sub='source48/input'if n==48 else f'transfers/dim{n}'
   run(f'{sub}/root_repair/check_tail.py');run(f'{sub}/root_repair/certify_annular.py')
  run('source48/check_exact.py');run('transfers/check_exact.py')
 print('QUICK EXACT CHECKS PASS; NOT FULL GEOMETRY REPLAY'if a.quick else f'COMPLETE SELECTED POWER6 LATTICE PROOF PASS {dims}',flush=True)
