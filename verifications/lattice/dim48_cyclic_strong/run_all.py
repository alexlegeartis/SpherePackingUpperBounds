import hashlib,subprocess,sys
from pathlib import Path
if not __debug__:raise RuntimeError('Assertions required')
p=Path(__file__).resolve().parent
assert hashlib.sha256((p/'certificate.json').read_bytes()).hexdigest()=='a40d47f8297b6475f8173ae8ba3f1c84356ecc36c43c18094b3cf31dc3c85df1'
for script,args in [('cyclic/meanvalue_angular.py',[p/'cyclic']),('check_all_short.py',[]),('check_narrow.py',[]),('replay.py',[p/'certificate.json']),('profiles.py',[p/'certificate.json']),('finish.py',[]),('verify_second.py',[p/'certificate.json'])]:
 subprocess.run([sys.executable,str(p/script),*map(str,args)],cwd=p,check=True)
print('FULL PORTABLE 48D PASS: density < 0.000000354687644',flush=True)
