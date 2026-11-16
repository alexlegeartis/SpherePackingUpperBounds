import hashlib,subprocess,sys
from pathlib import Path
if not __debug__:raise RuntimeError('Assertions required')
p=Path(__file__).resolve().parent
assert hashlib.sha256((p/'certificate.json').read_bytes()).hexdigest()=='5b6ae7366e992d6af08549ca0c29ea4f4e6a45bb30f5101bf76029d41341661f'
for script,args in [('cyclic/meanvalue_angular.py',[p/'cyclic']),('check_all_short.py',[]),('check_narrow.py',[]),('replay.py',[p/'certificate.json']),('profiles.py',[p/'certificate.json']),('finish.py',[]),('verify_second.py',[p/'certificate.json'])]:
 subprocess.run([sys.executable,str(p/script),*map(str,args)],cwd=p,check=True)
print('FULL PORTABLE 37D PASS: density < 0.000018573668945',flush=True)
