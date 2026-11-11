from pathlib import Path
import subprocess,sys,json,shutil
import numpy as np,mpmath as mp
from fractions import Fraction as F
from scipy.optimize import minimize_scalar
D=Path(__file__).resolve().parent;mp.mp.dps=120

def add_pair_cuts():
 d=json.loads((D/'unverified_export.json').read_text());cv=lambda s:mp.mpf(F(s).numerator)/F(s).denominator;p=list(map(cv,d['physical_polynomial']));M=cv(d['M_exact']);kap=mp.mpf(17)/8
 def value(x):
  X=mp.mpf(float(x));h=mp.exp(-kap*(X-6*mp.pi));return float(mp.polyval(list(reversed(p)),X)*mp.exp((kap-1)*X-kap*6*mp.pi)/(2-h)/M)
 grid=np.linspace(float(6*mp.pi),40,501);peaks=[]
 for l,u in zip(grid[:-1],grid[1:]):
  opt=minimize_scalar(lambda x:-value(x),bounds=(l,u),method='bounded',options={'xatol':1e-12})
  if -opt.fun>1:peaks.append((float(opt.x),-float(opt.fun)))
 extra=D/'radial_perturb46_degree47_extra.json';points=json.loads(extra.read_text());points['physical'] += [x for x,_ in peaks];extra.write_text(json.dumps(points,indent=2))
 return peaks
history=json.loads((D/'refinement_history.json').read_text());print('initial_pair_cuts',add_pair_cuts(),flush=True)
for i in range(8,13):
 search=subprocess.run([sys.executable,str(D/'radial_search.py'),'--dimension','46','--degree','47'],capture_output=True,text=True);(D/'history'/f'round{i}_search.log').write_text(search.stdout+search.stderr)
 if search.returncode:print('SEARCH_ERROR',i,search.stderr,flush=True);break
 candidate=D/'radial_perturb46_degree47.json';shutil.copy(candidate,D/'history'/f'round{i}_candidate.json')
 audit=subprocess.run([sys.executable,str(D/'audit_radial.py'),str(candidate)],capture_output=True,text=True);(D/'history'/f'round{i}_audit.log').write_text(audit.stdout+audit.stderr)
 row={'round':i,'density_candidate':json.loads(candidate.read_text())['density_candidate'],'radial_audit_returncode':audit.returncode}
 for kind in ['rejection','certificate']:
  p=D/f'radial_perturb46_degree47_{kind}.json'
  if p.exists():shutil.copy(p,D/'history'/f'round{i}_{kind}.json')
 if audit.returncode==0:row['radial_audit_pass']=True
 elif 'Forbidden sign roots' in audit.stderr:row['forbidden_root_counts']={k:len(v)for k,v in json.loads((D/'radial_perturb46_degree47_rejection.json').read_text())['global_sign_roots'].items()}
 elif "RuntimeError: coupling" in audit.stderr:row['added_pair_peaks']=add_pair_cuts()
 else:row['unexpected_failure']=audit.stderr
 history.append(row);print(row,flush=True);(D/'refinement_history.json').write_text(json.dumps(history,indent=2))
 if audit.returncode==0 or 'unexpected_failure'in row:break
