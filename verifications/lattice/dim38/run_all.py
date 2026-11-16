"""Full portable 38D replay. No network or optimization solver required."""
import hashlib,subprocess,sys
from pathlib import Path
if not __debug__:raise RuntimeError('Assertions required')
p=Path(__file__).resolve().parent
assert hashlib.sha256((p/'certificate.json').read_bytes()).hexdigest()=='42d13a9e04e2714f7768d55bcddbe45c6d3aa6c45e816190fd6cc5b38c89d96c'
for script,args in [('cyclic/primary_angular.py',[p/'cyclic']),('cyclic/independent_angular.py',[p/'cyclic']),('check_narrow38.py',[]),('check_all_short.py',[]),('replay38.py',[p/'certificate.json']),('replay38_profiles.py',[p/'certificate.json']),('verify_second.py',[p/'certificate.json'])]:
 subprocess.run([sys.executable,str(p/script),*map(str,args)],cwd=p,check=True)
print('FULL PORTABLE 38D LATTICE PASS: density < 0.000013234279694',flush=True)
