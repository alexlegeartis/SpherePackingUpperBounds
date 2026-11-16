"""Full portable 38D replay. No network or optimization solver required."""
import hashlib,subprocess,sys
from pathlib import Path
if not __debug__:raise RuntimeError('Assertions required')
p=Path(__file__).resolve().parent
assert hashlib.sha256((p/'certificate.json').read_bytes()).hexdigest()=='34888f6e673b2369bf69cd506e6fab0c3cf8612a03c3617ff85ca0d1678a9eda'
for script,args in [('cyclic/primary_angular.py',[p/'cyclic']),('cyclic/independent_angular.py',[p/'cyclic']),('check_narrow38.py',[]),('check_all_short.py',[]),('replay38.py',[p/'certificate.json']),('replay38_profiles.py',[p/'certificate.json']),('verify_second.py',[p/'certificate.json'])]:
 subprocess.run([sys.executable,str(p/script),*map(str,args)],cwd=p,check=True)
print('FULL PORTABLE 38D LATTICE PASS: density < 0.000013058679142',flush=True)
