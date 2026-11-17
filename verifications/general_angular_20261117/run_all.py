from pathlib import Path
import subprocess,sys,os,json,hashlib
P=Path(__file__).resolve().parent
if not __debug__:raise RuntimeError('Assertions required')
manifest=json.loads((P/'MANIFEST.json').read_text())
for f,h in manifest.items():assert hashlib.sha256((P/f).read_bytes()).hexdigest()==h,f
for n in [45,48]:
 for script in ['angular_structure_verify.py','angular_taylor_standalone.py','primary_radial.py','independent_radial.py']:
  subprocess.run([sys.executable,'-B',str(P/script),str(P/f'dim{n}')],check=True,env={**os.environ,'ANGULAR_PREC':'768'})
subprocess.run([sys.executable,'-B',str(P/'check_substitution.py')],check=True)
print('PASS complete reviewed unrestricted45/48 audit replay')
