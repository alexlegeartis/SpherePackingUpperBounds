from pathlib import Path
import subprocess,sys
assert __debug__
p=Path(__file__).resolve().parent
for component in ['angular','radial','domain']:
 subprocess.run([sys.executable,str(p/f'continuum_three_point/independent_cyclic40_{component}.py')],cwd=p,check=True)
subprocess.run([sys.executable,str(p/'verify_algebra.py')],check=True)
print('PASS independent47D restriction/Gram audit and full inherited source proofs')
