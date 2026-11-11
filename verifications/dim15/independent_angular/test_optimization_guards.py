"""Reject assertion-disabled executions before any proof checks or output writes."""
if not __debug__:
    raise RuntimeError('Assertions are required for the guard test.')
from pathlib import Path
import subprocess,sys,json
D=Path(__file__).resolve().parent;rows=[]
for script in ['verify_all.py','independent_check.py','primary_replay/verify_angular.py']:
 for flag in ['-O','-OO']:
  p=subprocess.run([sys.executable,flag,str(D/script)],capture_output=True,text=True)
  success=p.returncode!=0 and 'Assertions are required; do not run this verifier with -O or -OO.' in p.stderr
  if not success:raise RuntimeError(f'Guard did not reject {script} {flag}: {p.stdout} {p.stderr}')
  rows.append({'script':script,'flag':flag,'returncode':p.returncode,'rejected_as_expected':True,'stderr':p.stderr})
(D/'optimization_guard_results.json').write_text(json.dumps({'pass':True,'tests':rows},indent=2)+'\n');print('PASS: all three verifier scripts reject -O and -OO',flush=True)
