import hashlib,subprocess,sys
from pathlib import Path
if not __debug__:raise RuntimeError('Assertions required')
p=Path(__file__).resolve().parent
assert hashlib.sha256((p/'certificate.json').read_bytes()).hexdigest()=='bf8854ad2d867488b800cdb710d4db3085a324a012970854ede54866349bb55f'
for script,args in [('cyclic/primary_angular.py',[p/'cyclic']),('cyclic/independent_angular.py',[p/'cyclic'])]+[('replay.py',[p/'certificate.json']),('profiles.py',[p/'certificate.json']),('verify_second.py',[p/'certificate.json'])]:
 subprocess.run([sys.executable,str(p/script),*map(str,args)],cwd=p,check=True)
print('FULL PORTABLE 40D PASS: density < 0.000006869691783',flush=True)
