"""Hash-check a temporary copy and replay the complete selected proof chain."""
import argparse,hashlib,json,shutil,subprocess,sys,tempfile
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
if not __debug__:raise RuntimeError('Assertions required; do not use -O')
p=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--dimension',type=int,choices=range(36,49));ap.add_argument('--jobs',type=int,default=2);args=ap.parse_args();dims=[args.dimension]if args.dimension else list(range(36,49))
import gzip
stored=json.loads((p/'STORED_FILES.json').read_text());decoded={}
for name,digest in json.loads((p/'MANIFEST.json').read_text()).items():
 rec=stored.get(name);raw=(p/(rec['stored_path']if rec else name)).read_bytes()
 if rec:
  assert hashlib.sha256(raw).hexdigest()==rec['stored_sha256'],('stored hash mismatch',name)
  raw=gzip.decompress(raw);decoded[name]=raw
 assert hashlib.sha256(raw).hexdigest()==digest,('hash mismatch',name)
with tempfile.TemporaryDirectory(prefix='rank-two-release-')as td:
 root=Path(td)/'proof';shutil.copytree(p,root)
 for name,raw in decoded.items():(root/name).write_bytes(raw)
 def run(file,cwd=None,arguments=()):
  f=root/file;cwd=root if cwd is None else root/cwd
  result=subprocess.run([sys.executable,str(f),*map(str,arguments)],cwd=cwd,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  print(result.stdout,flush=True);assert result.returncode==0,('failed',file)
  if file.startswith('rank2_taylor'):assert 'DOMINANT GATE True' in result.stdout,('incomplete gate',file)
 def gate(n):
  k=37 if n<=37 else 42 if n<=42 else 43 if n==43 else 45 if n<=46 else 48
  return {37:'rank2_taylor37_audit100.py',42:'rank2_taylor42_recovery80.py',43:'rank2_taylor43_audit100.py',45:'rank2_taylor45_recovery80.py',48:'rank2_taylor_interval_audit60.py'}[k]
 tasks=sorted(set(gate(n)for n in dims))
 with ThreadPoolExecutor(max_workers=max(1,min(args.jobs,4)))as pool:list(pool.map(run,tasks))
 run('check_rank2_tails_exact.py')
 def dependency(n):
  if n==36:
   run('old_audit36_independent/verify_complete.py','old_audit36_independent')
   run('multiplier36_independent/verify.py','multiplier36_independent')
   run('multiplier36_independent/check_fourier_transform.py','multiplier36_independent');return
  folder='old_radial_audit'if n==48 else f'old_audit{n}'
  for name in ['check_all_short.py','check_narrow38.py'if n==38 else'check_narrow.py']:
   if(root/folder/name).exists():run(folder+'/'+name,folder)
  run(folder+('/replay38.py'if n==38 else'/replay.py'),folder,[root/folder/'certificate.json'])
  run(folder+('/replay38_profiles.py'if n==38 else'/profiles.py'),folder,[root/folder/'certificate.json'])
  if n in [46,48]:
   extra='multiplier46_independent'if n==46 else'combined48_independent'
   run(extra+'/verify.py',extra);run(extra+'/check_fourier_transform.py',extra)
   if n==48:run(extra+'/check_sos.py',extra)
 with ThreadPoolExecutor(max_workers=max(1,min(args.jobs,4)))as pool:list(pool.map(dependency,dims))
 if any(n in dims for n in [46,47]):
  run('old_radial_audit/check_all_short.py','old_radial_audit');run('old_radial_audit/check_narrow.py','old_radial_audit')
 run('check_all_transferred_bounds.py');run('check_transfer37.py')
 rows={r['dimension']:r for r in json.loads((root/'VERIFIED_RESULTS.json').read_text())}
 registered=json.loads((root/'REGISTERED_RESULTS.json').read_text())
 for r in registered:
  n=r['dimension'];safe=Fraction(r['upper'])
  if n in [46,48]:
   extra='multiplier46_independent'if n==46 else'combined48_independent'
   if n in dims:assert Fraction(json.loads((root/extra/'certificate.json').read_text())['density'])<safe
  elif n>=38:assert Fraction(rows[n]['verified_upper'])<=safe
  elif n==37:assert Fraction(json.loads((root/'PROSPECTIVE_RESULTS_37.json').read_text())[0]['safe'])<=safe
  elif n in dims:assert Fraction(json.loads((root/'multiplier36_independent/certificate.json').read_text())['density'])<safe
 print('COMPLETE SELECTED LATTICE PROOF PASS',dims,flush=True)
