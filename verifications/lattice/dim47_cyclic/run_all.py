import hashlib,subprocess,sys
from pathlib import Path
if not __debug__:raise RuntimeError('Assertions required')
p=Path(__file__).resolve().parent
assert hashlib.sha256((p/'certificate.json').read_bytes()).hexdigest()=='e861afbdbe240be629b42d3b96c85c1db19e327e6e55e1241bad4c487121df76'
for script,args in [('cyclic/meanvalue_angular.py',[p/'cyclic']),('check_all_short.py',[]),('check_narrow.py',[]),('replay.py',[p/'certificate.json']),('profiles.py',[p/'certificate.json']),('finish.py',[]),('verify_second.py',[p/'certificate.json'])]:
 subprocess.run([sys.executable,str(p/script),*map(str,args)],cwd=p,check=True)
print('FULL PORTABLE 47D PASS: density < 0.000000527909479',flush=True)
