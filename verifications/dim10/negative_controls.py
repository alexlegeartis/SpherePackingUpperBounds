"""Mutation tests: failures must remain failures, without changing proof inputs."""
from pathlib import Path
import tempfile,shutil,subprocess,sys,json
D=Path(__file__).resolve().parent
cases=[('negative_angular','angular.json','gegenbauer_coefficients',lambda v:['-1']+v[1:],'verify_geometry_independent.py'),('negative_profile','fit.json','parameters',lambda v:[-1]+v[1:],'verify_geometry_independent.py'),('wrong_fourier','certificate.json','fourier',lambda v:['2']+v[1:],'verify_radial_independent.py'),('zero_C','certificate.json','C_exact',lambda v:'0','verify_radial_independent.py'),('missing_cover_leaf','partition.json','accepted',lambda v:v[1:],'verify_geometry_independent.py')]
for name,file,key,change,script in cases:
 with tempfile.TemporaryDirectory() as tmp:
  target=Path(tmp)/'audit';shutil.copytree(D/'independent',target)
  obj=json.loads((target/file).read_text());obj[key]=change(obj[key]);(target/file).write_text(json.dumps(obj))
  r=subprocess.run([sys.executable,str(target/script)],capture_output=True,text=True,timeout=120)
  assert r.returncode!=0,(name,'mutation accepted')
 print('PASS rejected',name,flush=True)
for script in ['verify_all.py','geometry/verify_geometry.py','radial/verify.py','independent/verify_geometry_independent.py','independent/verify_radial_independent.py']:
 r=subprocess.run([sys.executable,'-O',str(D/script)],capture_output=True,text=True,timeout=10)
 assert r.returncode!=0,('optimized_python',script)
print('PASS: 5 mathematical mutation controls and 5 optimized-Python rejection controls')
