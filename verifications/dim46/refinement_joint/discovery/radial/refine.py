from pathlib import Path
import subprocess,sys,json,shutil
D=Path(__file__).resolve().parent
history=[]
for i in range(1,9):
 search=subprocess.run([sys.executable,str(D/'radial_search.py'),'--dimension','46','--degree','31'],capture_output=True,text=True)
 (D/'history'/f'round{i}_search.log').write_text(search.stdout+search.stderr)
 if search.returncode:print('SEARCH_ERROR',i,search.stderr,flush=True);break
 candidate=D/'radial_perturb46_degree31.json';shutil.copy(candidate,D/'history'/f'round{i}_candidate.json')
 audit=subprocess.run([sys.executable,str(D/'audit_radial.py'),str(candidate)],capture_output=True,text=True)
 (D/'history'/f'round{i}_audit.log').write_text(audit.stdout+audit.stderr)
 row={'round':i,'density_candidate':json.loads(candidate.read_text())['density_candidate'],'radial_audit_returncode':audit.returncode}
 for kind in ['rejection','certificate']:
  p=D/f'radial_perturb46_degree31_{kind}.json'
  if p.exists():shutil.copy(p,D/'history'/f'round{i}_{kind}.json')
 if audit.returncode==0:row['radial_audit_pass']=True
 else:
  if 'Forbidden sign roots' not in audit.stderr:row['unexpected_failure']=audit.stderr
  else:row['forbidden_root_counts']={k:len(v)for k,v in json.loads((D/'radial_perturb46_degree31_rejection.json').read_text())['global_sign_roots'].items()}
 history.append(row);print(row,flush=True);(D/'refinement_history.json').write_text(json.dumps(history,indent=2))
 if audit.returncode==0 or 'unexpected_failure'in row:break
