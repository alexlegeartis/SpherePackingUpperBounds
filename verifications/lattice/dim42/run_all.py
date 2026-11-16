import hashlib,subprocess,sys
from pathlib import Path
if not __debug__:raise RuntimeError('Assertions required')
p=Path(__file__).resolve().parent
assert hashlib.sha256((p/'certificate.json').read_bytes()).hexdigest()=='bc1334726e33164aa41da778590890461e782db2cd50c34ebc63aae6e1a25c7d'
for script,args in [('cyclic/meanvalue_angular.py',[p/'cyclic']),('replay.py',[p/'certificate.json']),('profiles.py',[p/'certificate.json']),('verify_second.py',[p/'certificate.json'])]:
 subprocess.run([sys.executable,str(p/script),*map(str,args)],cwd=p,check=True)
print('FULL PORTABLE 42D PASS: density < 0.000003571509824',flush=True)
