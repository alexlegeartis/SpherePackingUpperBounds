"""Re-export frozen discovery checkpoints in isolation and compare exact bytes."""
from pathlib import Path
import hashlib,json,os,shutil,subprocess,sys,tempfile,time
if not __debug__:raise RuntimeError('Assertions required')
D=Path(__file__).resolve().parent
manifest=json.loads((D/'INPUT_MANIFEST.json').read_text())
expected=json.loads((D/'EXPECTED_EXPORTS.json').read_text())
for name,record in manifest.items():
 assert hashlib.sha256((D/name).read_bytes()).hexdigest()==record['sha256']
t0=time.monotonic();runs=[]
with tempfile.TemporaryDirectory(prefix='packing22_23_export_')as tmp:
 W=Path(tmp)
 for name in ['inputs','templates']:shutil.copytree(D/name,W/name)
 for name in ['export_angular_n.py']:shutil.copy2(D/name,W/name)
 for n in [22,23]:
  tag='2p75'if n==22 else'2p875';candidate=f'radial_scaled{n}_degree31_r{tag}.json'
  for name in [f'angular{n}_a40_scaled_candidate.json',candidate]:shutil.copy2(D/name,W/name)
  commands=[[sys.executable,'export_angular_n.py',str(n)],
   [sys.executable,f'dim{n}/independent_cyclic40_angular.py'],
   [sys.executable,f'dim{n}/independent_cyclic40_domain.py'],
   [sys.executable,f'dim{n}/audit_radial.py',candidate]]
  for cmd in commands:
   result=subprocess.run(cmd,cwd=W,text=True,capture_output=True,env={**os.environ,'OPENBLAS_NUM_THREADS':'1'},timeout=180)
   print(result.stdout,end='',flush=True)
   if result.returncode:print(result.stderr,flush=True);raise RuntimeError(f'Failed: {cmd}')
   runs.append({'command':cmd[1:],'pass':True})
 for name,record in expected.items():
  actual=hashlib.sha256((W/name).read_bytes()).hexdigest();assert actual==record['sha256'],(name,actual,record['sha256'])
  print('BYTE_IDENTICAL',name,flush=True)
report={'pass':True,'scope':'Frozen checkpoint export and full angular/radial replay; no optimizer rerun.','all_four_exports_byte_identical':True,'input_hashes_match':True,'runs':runs,'seconds':time.monotonic()-t0}
(D/'frozen_export_replay.json').write_text(json.dumps(report,indent=2)+'\n');print(report,flush=True)
