import hashlib,subprocess,sys
from pathlib import Path
if not __debug__:raise RuntimeError('Assertions required')
p=Path(__file__).resolve().parent
assert hashlib.sha256((p/'certificate.json').read_bytes()).hexdigest()=='1e2c449de10236a22adbab00fef0ffcac4bf82e5bfb9128e16b8336b34060ea0'
for script,args in [('cyclic/meanvalue_angular.py',[p/'cyclic']),('replay.py',[p/'certificate.json']),('profiles.py',[p/'certificate.json']),('verify_second.py',[p/'certificate.json'])]:
 subprocess.run([sys.executable,str(p/script),*map(str,args)],cwd=p,check=True)
print('FULL PORTABLE 42D PASS: density < 0.000003536520866',flush=True)
