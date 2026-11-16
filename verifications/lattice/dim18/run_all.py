"""Portable18D complete proof replay. Run without -O; Python3+SymPy."""
from pathlib import Path
import subprocess,sys
if not __debug__:raise RuntimeError('Do not run with -O')
p=Path(__file__).resolve().parent
jobs=[[p.parent/'prior'/'taylor.py',p/'prior_certificate.json'],[p/'verify_determinant_thresholds.py'],[p/'verify_small_torsion.py'],[p/'verify_remaining_torsion.py'],[p/'verify_certificate_independent.py']]
for args in jobs:
 print('RUN',args[0].name,flush=True)
 subprocess.run([sys.executable,*map(str,args)],cwd=p,check=True)
print('COMPLETE18D PRIOR, TRACE THEOREM AND DENSITY CERTIFICATE PASS',flush=True)
