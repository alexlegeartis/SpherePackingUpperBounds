import hashlib,subprocess,sys
from pathlib import Path
if not __debug__:raise RuntimeError('Assertions required')
p=Path(__file__).resolve().parent
assert hashlib.sha256((p/'certificate.json').read_bytes()).hexdigest()=='6c439b6ee1d2afae624f21587c8ad79cc715985b741b4c7b9c1234e67cb86239'
for script,args in [('cyclic/meanvalue_angular.py',[p/'cyclic']),('replay.py',[p/'certificate.json']),('profiles.py',[p/'certificate.json']),('verify_second.py',[p/'certificate.json'])]:
 subprocess.run([sys.executable,str(p/script),*map(str,args)],cwd=p,check=True)
print('FULL PORTABLE 43D PASS: density < 0.000002566728865',flush=True)
