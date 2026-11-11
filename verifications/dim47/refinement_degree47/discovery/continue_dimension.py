"""Up to six further root-guided repairs, preserving prior checkpoints."""
from pathlib import Path
import subprocess,sys,json,shutil,argparse,time
D=Path(__file__).resolve().parent;ap=argparse.ArgumentParser();ap.add_argument('dimension',type=int,choices=[46,47]);a=ap.parse_args();n=a.dimension;history=json.loads((D/f'history{n}.json').read_text());start=time.time();base=json.loads((D/f'baseline{n}_degree31.json').read_text());H=D/f'history{n}'
for i in range(5,11):
 row={'round':i,'repair_number':i,'degree':47,'dimension':n}
 try:search=subprocess.run([sys.executable,str(D/'radial_search.py'),'--dimension',str(n),'--degree','47','--control-safe-margin'],capture_output=True,text=True,timeout=180)
 except subprocess.TimeoutExpired as exc:
  row.update({'search_timeout':True,'error':str(exc)});history.append(row);break
 (H/f'round{i}_search.log').write_text(search.stdout+search.stderr)
 if search.returncode:
  row.update({'search_returncode':search.returncode,'search_failed':True,'error':search.stderr or search.stdout});history.append(row);print(row,flush=True);break
 candidate=D/f'radial_perturb{n}_degree47.json';shutil.copy2(candidate,H/f'round{i}_candidate.json');can=json.loads(candidate.read_text());row.update({'density_candidate':can['density_candidate'],'relative_numerical_gain_vs_degree31':1-can['density_candidate']/base['density_candidate']})
 try:audit=subprocess.run([sys.executable,str(D/'audit_radial.py'),str(candidate)],capture_output=True,text=True,timeout=180)
 except subprocess.TimeoutExpired as exc:
  row.update({'audit_timeout':True,'error':str(exc)});history.append(row);print(row,flush=True);break
 (H/f'round{i}_audit.log').write_text(audit.stdout+audit.stderr);row['radial_audit_returncode']=audit.returncode
 exact=D/f'radial_perturb{n}_degree47_radial_exact.json'
 if exact.exists():shutil.copy2(exact,H/f'round{i}_radial_exact.json')
 if audit.returncode==0:
  cert=D/f'radial_perturb{n}_degree47_certificate.json';shutil.copy2(cert,H/f'round{i}_certificate.json');row['radial_audit_pass']=True;row['density_interval']=json.loads(cert.read_text())['density_interval'];row['complete_package_audit_pass']=False
 elif 'Forbidden sign roots'in audit.stderr:
  reject=D/f'radial_perturb{n}_degree47_rejection.json';shutil.copy2(reject,H/f'round{i}_rejection.json');row['forbidden_root_counts']={k:len(v)for k,v in json.loads(reject.read_text())['global_sign_roots'].items()}
 else:row['unexpected_failure']=audit.stderr
 history.append(row);(D/f'history{n}.json').write_text(json.dumps(history,indent=2)+'\n');print(row,flush=True)
 if audit.returncode==0 or 'unexpected_failure'in row:break
 if row['relative_numerical_gain_vs_degree31']<=1e-6:row['stop_reason']='Numerical improvement no greater than one part per million; no more repairs.';break
(D/f'history{n}.json').write_text(json.dumps(history,indent=2)+'\n');print({'dimension':n,'seconds':time.time()-start,'completed_solve_attempts':len(history),'complete_package_audit_pass':False},flush=True)
