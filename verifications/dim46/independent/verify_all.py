from pathlib import Path
import subprocess,sys,hashlib
assert __debug__
p=Path(__file__).resolve().parent
assert (p/'angular.json').read_bytes()==(p/'continuum_three_point/cyclic_angular40_exact.json').read_bytes()
for component in ['angular','domain']:
 subprocess.run([sys.executable,str(p/f'continuum_three_point/independent_cyclic40_{component}.py')],cwd=p,check=True)
subprocess.run([sys.executable,str(p/'verify.py')],check=True)
print('PASS independent46D radial and inherited angular/cyclic proof chain')
