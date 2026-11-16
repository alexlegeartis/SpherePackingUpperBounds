from pathlib import Path
import subprocess,sys
if not __debug__:raise RuntimeError('Do not run with python -O')
p=Path(__file__).resolve().parent
for args in [[p/'run_theorem.py'],[p/'verify_certificate_independent.py']]:
 print('RUN',args[0].name,flush=True)
 subprocess.run([sys.executable,*map(str,args)],cwd=p,check=True)
print('COMPLETE 14D LATTICE CERTIFICATE PASS')
