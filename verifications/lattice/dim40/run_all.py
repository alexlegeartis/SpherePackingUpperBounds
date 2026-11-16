"""Full portable 40D lattice replay; no optimizer or network."""
import hashlib,subprocess,sys
from pathlib import Path
if not __debug__:raise RuntimeError('Assertions required')
p=Path(__file__).resolve().parent
assert hashlib.sha256((p/'certificate.json').read_bytes()).hexdigest()=='01b7ab491970e717c9dec32f13f42fdf326bb26826518d509dae22f34fb4ff28'
for script,args in [('cyclic/primary_angular.py',[p/'cyclic']),('cyclic/independent_angular.py',[p/'cyclic']),('replay40.py',[p/'certificate.json']),('replay40_profiles.py',[p/'certificate.json']),('verify_second.py',[p/'certificate.json'])]:
 subprocess.run([sys.executable,str(p/script),*map(str,args)],cwd=p,check=True)
print('FULL PORTABLE 40D LATTICE PASS: density < 0.000006887626280',flush=True)
