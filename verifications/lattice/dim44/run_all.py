import hashlib,subprocess,sys
from pathlib import Path
if not __debug__:raise RuntimeError('Assertions required')
p=Path(__file__).resolve().parent
assert hashlib.sha256((p/'certificate.json').read_bytes()).hexdigest()=='2b373f4af6df2b750836328262f79af78f2f72ef7ccfd00e659913c3b3ad923d'
for script,args in [('cyclic/meanvalue_angular.py',[p/'cyclic']),('replay.py',[p/'certificate.json']),('profiles.py',[p/'certificate.json']),('verify_second.py',[p/'certificate.json'])]:
 subprocess.run([sys.executable,str(p/script),*map(str,args)],cwd=p,check=True)
print('FULL PORTABLE 44D PASS: density < 0.000001783944991',flush=True)
