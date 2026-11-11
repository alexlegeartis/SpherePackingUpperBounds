"""Bounded fresh discovery in a new output folder; no record claim."""
from pathlib import Path
import argparse,os,shutil,subprocess,sys
if not __debug__:raise RuntimeError("Assertions required")
D=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--dimension',type=int,choices=[22,23],required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--repair-rounds',type=int,default=10);args=ap.parse_args()
assert 1<=args.repair_rounds<=20
W=args.output.resolve()
if W.exists():raise RuntimeError('Output must not already exist; protect frozen checkpoints')
W.mkdir(parents=True)
for name in ['inputs','templates']:shutil.copytree(D/name,W/name)
for name in ['search_angular.py','export_angular_n.py','search_robust.py']:shutil.copy2(D/name,W/name)
n=args.dimension;r2='2.75'if n==22 else'2.875';tag=r2.replace('.','p')
def run(arguments,label,timeout=180):
 p=subprocess.run([sys.executable,*arguments],cwd=W,env={**os.environ,'OPENBLAS_NUM_THREADS':'1'},text=True,capture_output=True,timeout=timeout)
 (W/f'{label}.log').write_text(p.stdout+p.stderr);print(label,p.returncode,flush=True);return p.returncode
assert run(['search_angular.py','--dimension',str(n),'--degree','19','--alpha','40','--rounds','8'],'search_angular')==0
assert run(['export_angular_n.py',str(n)],'export_angular')==0
for part in ['angular','domain']:assert run([f'dim{n}/independent_cyclic40_{part}.py'],f'verify_{part}')==0
for iteration in range(1,args.repair_rounds+1):
 assert run(['search_robust.py','--dimension',str(n),'--degree','31','--r2',r2],f'search_radial_{iteration}')==0
 if run([f'dim{n}/audit_radial.py',f'radial_scaled{n}_degree31_r{tag}.json'],f'audit_radial_{iteration}')==0:
  print('DISCOVERY_EXACT_CHECKS_PASS',W,flush=True);break
else:raise RuntimeError('Bounded repair budget exhausted; no certified result')
