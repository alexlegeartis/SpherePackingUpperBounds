"""Complete independently reviewed 7D frozen count certificate audit, no optimizer."""
from pathlib import Path
import subprocess,sys
if not __debug__:raise RuntimeError('Do not disable assertions')
p=Path(__file__).resolve().parent
for script in ['verify_algebra_guards.py','verify_pair.py','certify_far_bernstein.py','certify_cyclic_bernstein.py']:
 args=[sys.executable,str(p/script)]
 if script.startswith('certify_'):args.append(str(p/'witness_psd_exact.json'))
 subprocess.run(args,cwd=p,check=True)
print('PASS complete independently reviewed frozen seven-dimensional count certificate')
