from pathlib import Path
import subprocess,sys
if not __debug__:raise RuntimeError('Do not run with python -O')
p=Path(__file__).resolve().parent
for args in [[p.parent/'prior'/'taylor.py',p/'prior_certificate.json'],[p/'verify_determinant.py'],[p/'verify_arithmetic.py'],[p/'verify_certificate_independent.py']]:
 print('RUN',args[0].name,flush=True)
 subprocess.run([sys.executable,*map(str,args)],cwd=p,check=True)
print('COMPLETE 15D LATTICE CERTIFICATE PASS')
