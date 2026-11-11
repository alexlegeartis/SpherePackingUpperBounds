"""Bounded root-guided repair. All rejected checkpoints are retained."""
from pathlib import Path
import subprocess,json,sys,shutil,time
D=Path(__file__).resolve().parent;tag='radial_cut12';hist=[];H=D/'history_cut12';H.mkdir(exist_ok=True)
for rnd in range(1,13):
 j=json.loads((D/(tag+'.json')).read_text());cutsfile=D/(tag+'_cuts.json');cuts=json.loads(cutsfile.read_text()) if cutsfile.exists() else {'physical':[],'fourier':[],'pair':[]}
 reject=D/(tag+'_rejection.json')
 if reject.exists():
  for key,roots in json.loads(reject.read_text())['global_sign_roots'].items():
   xx=sorted(x[0]for x in roots);cuts.setdefault(key,[]).extend(xx+[(x+y)/2 for x,y in zip(xx,xx[1:])])
 pair=D/(tag+'_pair_rejection.json')
 if pair.exists():cuts.setdefault('pair',[]).append(json.loads(pair.read_text())['X'])
 cutsfile.write_text(json.dumps(cuts,indent=2)+'\n')
 for suff in ['.json','_rejection.json','_pair_rejection.json','_exact.json']:
  f=D/(tag+suff)
  if f.exists():shutil.copy(f,H/f'round{rnd-1}{suff}')
 if pair.exists():pair.unlink()
 cmd=[sys.executable,str(D/'search_radial.py'),'--dimension','12','--degree','63','--s','12/8','--a','10','--b','16','--tag',tag,'--margin','1e-8','--cut','12']
 with (H/f'search{rnd}.log').open('w')as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=150)
 assert r.returncode==0
 new=json.loads((D/(tag+'.json')).read_text());print('Round',rnd,'candidate',new['density_candidate'],flush=True)
 with (H/f'audit{rnd}.log').open('w')as f:r=subprocess.run([sys.executable,str(D/'audit_radial.py'),str(D/(tag+'.json'))],stdout=f,stderr=subprocess.STDOUT,timeout=180)
 hist.append({'round':rnd,'density_candidate':new['density_candidate'],'radial_audit_pass':r.returncode==0});(D/'refinement_cut12_history.json').write_text(json.dumps(hist,indent=2)+'\n')
 if r.returncode==0:print('RADIAL PASS',flush=True);break
