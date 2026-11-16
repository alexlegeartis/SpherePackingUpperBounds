"""Full portable 39D lattice replay; no optimizer or network."""
import hashlib,subprocess,sys
from pathlib import Path
if not __debug__:raise RuntimeError('Assertions required')
p=Path(__file__).resolve().parent
assert hashlib.sha256((p/'certificate.json').read_bytes()).hexdigest()=='32915e94df222f092af74294917a4896df8a3f85d6c20db17e1ba085e368bc31'
for script,args in [('cyclic/primary_angular.py',[p/'cyclic']),('cyclic/independent_angular.py',[p/'cyclic']),('replay39.py',[p/'certificate.json']),('replay39_profiles.py',[p/'certificate.json']),('verify_second.py',[p/'certificate.json'])]:
 subprocess.run([sys.executable,str(p/script),*map(str,args)],cwd=p,check=True)
print('FULL PORTABLE 39D LATTICE PASS: density < 0.000009086501484',flush=True)
