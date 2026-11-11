"""Bounded root-guided repair. All rejected checkpoints are retained."""
from pathlib import Path
import subprocess,json,sys,shutil,time,math
from flint import fmpq,arb,arb_poly,ctx
from scipy.optimize import minimize_scalar
ctx.prec=384
if not __debug__:raise RuntimeError('Do not use Python -O')
D=Path(__file__).resolve().parent;tag='radial_cut12';hist=json.loads((D/'refinement_cut12_history.json').read_text());H=D/'history_cut12';H.mkdir(exist_ok=True)
for rnd in range(9,13):
 j=json.loads((D/(tag+'.json')).read_text());cutsfile=D/(tag+'_cuts.json');cuts=json.loads(cutsfile.read_text()) if cutsfile.exists() else {'physical':[],'fourier':[],'pair':[]}
 reject=D/(tag+'_rejection.json')
 if reject.exists():
  for key,roots in json.loads(reject.read_text())['global_sign_roots'].items():
   xx=sorted(x[0]for x in roots);cuts.setdefault(key,[]).extend(xx+[(x+y)/2 for x,y in zip(xx,xx[1:])])
 pair=D/(tag+'_pair_rejection.json')
 if pair.exists():
  bad=json.loads(pair.read_text())['X'];cuts.setdefault('pair',[]).append(bad)
  exact=json.loads((D/(tag+'_exact.json')).read_text());AP=arb_poly([arb(fmpq(x))for x in exact['physical']]);M=arb(fmpq(exact['M_exact']));ss=arb(fmpq(exact['s']));aa=arb(fmpq(exact['a_exact']));bb=arb(fmpq(exact['b_exact']));cut=float(exact['cut_X'])
  def gap(x):
   xx=arb(x);z=xx/(arb.pi()*ss)-1;hh=(-aa*z-bb*z*z).exp();return AP(xx)*(-xx).exp()-M*(2*hh-hh*hh)
  op=minimize_scalar(lambda x:-float(gap(x).mid()),bounds=(max(math.pi*float(ss.mid()),bad-.05),min(cut,bad+.05)),method='bounded',options={'xatol':1e-13,'maxiter':100})
  peak=float(op.x);newpoints=[peak+d for d in [0,-.00025,.00025,-.0005,.0005,-.001,.001,-.002,.002,-.005,.005]];cuts['pair'].extend(newpoints)
  (H/f'pair_peak{rnd}.json').write_text(json.dumps({'previous_first_positive':bad,'peak_X':peak,'peak_residual_interval':str(gap(peak)),'new_pair_cuts':newpoints},indent=2)+'\n')
 cutsfile.write_text(json.dumps(cuts,indent=2)+'\n')
 for suff in ['.json','_rejection.json','_pair_rejection.json','_exact.json']:
  f=D/(tag+suff)
  if f.exists():shutil.copy(f,H/f'round{rnd-1}{suff}')
 if pair.exists():pair.unlink()
 cmd=[sys.executable,str(D/'search_radial.py'),'--dimension','14','--degree','63','--s','7/4','--a','14','--b','8','--tag',tag,'--margin','1e-8','--cut','12']
 with (H/f'search{rnd}.log').open('w')as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=150)
 assert r.returncode==0
 new=json.loads((D/(tag+'.json')).read_text());print('Round',rnd,'candidate',new['density_candidate'],flush=True)
 assert new['density_candidate']<0.04501031,'Stopped: numerical lead lost competitiveness'
 with (H/f'audit{rnd}.log').open('w')as f:r=subprocess.run([sys.executable,str(D/'audit_radial.py'),str(D/(tag+'.json'))],stdout=f,stderr=subprocess.STDOUT,timeout=180)
 hist.append({'round':rnd,'density_candidate':new['density_candidate'],'radial_audit_pass':r.returncode==0});(D/'refinement_cut12_history.json').write_text(json.dumps(hist,indent=2)+'\n')
 if r.returncode==0:print('RADIAL PASS',flush=True);break
