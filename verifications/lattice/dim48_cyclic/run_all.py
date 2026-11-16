import hashlib,subprocess,sys
from pathlib import Path
if not __debug__:raise RuntimeError('Assertions required')
p=Path(__file__).resolve().parent
assert hashlib.sha256((p/'certificate.json').read_bytes()).hexdigest()=='a03938f8d309fcc2c237bd0df6342f45c44fced0765a1b453a66f85f5f6bad6f'
for script,args in [('cyclic/meanvalue_angular.py',[p/'cyclic']),('check_all_short.py',[]),('check_narrow.py',[]),('replay.py',[p/'certificate.json']),('profiles.py',[p/'certificate.json']),('finish.py',[]),('verify_second.py',[p/'certificate.json'])]:
 subprocess.run([sys.executable,str(p/script),*map(str,args)],cwd=p,check=True)
print('FULL PORTABLE 48D PASS: density < 0.000000358368032',flush=True)
