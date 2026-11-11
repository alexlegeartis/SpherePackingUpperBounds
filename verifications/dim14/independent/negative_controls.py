"""Tampered mathematical inputs and disabled assertions must fail closed."""
from pathlib import Path
from fractions import Fraction as F
import tempfile,shutil,json,subprocess,sys,hashlib
if not __debug__:raise RuntimeError('Do not use -O')
D=Path(__file__).resolve().parent;results=[]
for kind in ['disabled_angular_assertions','disabled_radial_assertions','physical_coefficient','normalization','origin_amplitude','negative_harmonic']:
 with tempfile.TemporaryDirectory() as tmp:
  t=Path(tmp);shutil.copytree(D/'inputs',t/'inputs')
  for script in ['audit_angular.py','audit_radial.py']:shutil.copy2(D/script,t/script)
  script='audit_angular.py' if kind in ['disabled_angular_assertions','negative_harmonic']else 'audit_radial.py'
  args=[sys.executable]
  if kind.startswith('disabled'):args+=['-O']
  elif kind=='negative_harmonic':
   p=t/'inputs/angular_exact.json';old=hashlib.sha256(p.read_bytes()).hexdigest();j=json.loads(p.read_text());j['gegenbauer_coefficients'][1]='-1';p.write_text(json.dumps(j));new=hashlib.sha256(p.read_bytes()).hexdigest()
   # Adjust the test-copy hash pin, so rejection tests mathematics, not integrity.
   c=t/script;c.write_text(c.read_text().replace(old,new))
  else:
   p=t/'inputs/certificate.json';j=json.loads(p.read_text())
   if kind=='physical_coefficient':j['physical'][0]=str(F(j['physical'][0])+1)
   elif kind=='normalization':j['normalization_exact']='-1'
   elif kind=='origin_amplitude':j['M_exact']=str(F(j['M_exact'])*2)
   p.write_text(json.dumps(j))
  run=subprocess.run(args+[str(t/script)],capture_output=True,text=True,timeout=120)
  assert run.returncode!=0,(kind,'UNEXPECTED ACCEPTANCE')
  results.append({'control':kind,'rejected':True,'returncode':run.returncode,'last_error':run.stderr.strip().splitlines()[-1]})
(D/'negative_controls_result.json').write_text(json.dumps({'pass':True,'controls':results},indent=2)+'\n')
print('NEGATIVE_CONTROLS_PASS',len(results))
