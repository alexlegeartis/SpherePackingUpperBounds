import hashlib,subprocess,sys
from pathlib import Path
if not __debug__:raise RuntimeError('Assertions required')
p=Path(__file__).resolve().parent
assert hashlib.sha256((p/'certificate.json').read_bytes()).hexdigest()=='8e9cd92e9595e8dfe28e1e19daa468829547119883aa4f7294c5cfa3bb0feb29'
for script,args in [('cyclic/meanvalue_angular.py',[p/'cyclic']),('check_all_short.py',[]),('check_narrow.py',[]),('replay.py',[p/'certificate.json']),('profiles.py',[p/'certificate.json']),('finish.py',[]),('verify_second.py',[p/'certificate.json'])]:
 subprocess.run([sys.executable,str(p/script),*map(str,args)],cwd=p,check=True)
print('FULL PORTABLE 45D PASS: density < 0.000001185930816',flush=True)
