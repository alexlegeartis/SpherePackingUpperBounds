if not __debug__:
 raise RuntimeError('Do not run with python -O')
import subprocess,sys
from pathlib import Path
p=Path(__file__).resolve().parent
jobs=[['verify_prior.py'],['verify_nearframe.py']]+[['verify_metric_families.py',mode] for mode in ['five','seven','nine','ternary']]+[['verify_prime_support_trees.py'],['verify_strong_constants.py'],['verify.py']]
for job in jobs:
 subprocess.run([sys.executable,str(p/job[0]),*job[1:]],check=True,cwd=p)
print('ALL CHECKS PASS: 16D lattice density <=0.018655117428')
