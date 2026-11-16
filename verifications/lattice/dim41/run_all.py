import hashlib,subprocess,sys
from pathlib import Path
if not __debug__:raise RuntimeError('Assertions required')
p=Path(__file__).resolve().parent
assert hashlib.sha256((p/'certificate.json').read_bytes()).hexdigest()=='087f3bae78580cfc01aa26c6e3ff2d0e31038bf095b6806121172b9ef35c1773'
for script,args in [('cyclic/primary_angular.py',[p/'cyclic']),('cyclic/meanvalue_angular.py',[p/'cyclic'])]+[('replay.py',[p/'certificate.json']),('profiles.py',[p/'certificate.json']),('verify_second.py',[p/'certificate.json'])]:
 subprocess.run([sys.executable,str(p/script),*map(str,args)],cwd=p,check=True)
print('FULL PORTABLE 41D PASS: density < 0.000004531256243',flush=True)
