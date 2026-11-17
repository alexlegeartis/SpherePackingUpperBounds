"""Hash-checked portable exact/interval replay in temporary storage."""
import argparse,json,hashlib,sys,subprocess,tempfile,gzip
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
if not __debug__:raise RuntimeError('Do not use python -O')
D=Path(__file__).resolve().parent;available=sorted(map(int,json.loads((D/'VERIFIED_RESULTS.json').read_text())['bounds']))
ap=argparse.ArgumentParser();ap.add_argument('--quick',action='store_true',help='Exact algebra/cap/tail only, NOT full geometry replay');ap.add_argument('--dimension',type=int,choices=available);ap.add_argument('--jobs',type=int,choices=[1,2],default=2);a=ap.parse_args();dims=[a.dimension]if a.dimension else available
manifest=json.loads((D/'MANIFEST.json').read_text());stored=json.loads((D/'STORED_FILES.json').read_text())
with tempfile.TemporaryDirectory(prefix='adapted-lattice-proof-')as td:
 P=Path(td)/'proof';P.mkdir()
 for f,h in manifest.items():
  entry=stored.get(f,{});sf=entry.get('stored_path',f)
  for rel in [f,sf]:
   rp=Path(rel)
   if rp.is_absolute()or'..'in rp.parts:raise ValueError('Unsafe manifest path')
  data=(D/sf).read_bytes()
  if entry:
   assert hashlib.sha256(data).hexdigest()==entry['stored_sha256'],('stored hash mismatch',f)
   data=gzip.decompress(data)
  assert hashlib.sha256(data).hexdigest()==h,('raw hash mismatch',f)
  target=P/f;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
 def run(f,*args):
  r=subprocess.run([sys.executable,str(P/f),*map(str,args)],cwd=P,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);print(r.stdout,flush=True);assert r.returncode==0,('failed',f)
 def verify(n):
  prefix=f'dim{n}';run(prefix+'/check_exact.py');run(prefix+'/cap_independent.py');run(prefix+'/check_tail_exact.py')
  if not a.quick:
   run(prefix+'/replay_trace.py','main');run(prefix+'/input/root_repair/check_tail.py');run(prefix+'/input/root_repair/certify_annular.py');run(prefix+'/check_exact.py')
  print('DIMENSION',n,'QUICK ONLY'if a.quick else'FULL REPLAY PASS',flush=True)
 with ThreadPoolExecutor(max_workers=a.jobs)as pool:list(pool.map(verify,dims))
 print('QUICK EXACT CHECKS PASS; NOT FULL GEOMETRY REPLAY'if a.quick else f'COMPLETE ADAPTED LATTICE PROOF PASS {dims}',flush=True)
