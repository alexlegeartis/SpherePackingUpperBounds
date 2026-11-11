from pathlib import Path
import subprocess,sys
if not __debug__: raise RuntimeError('Assertions must be enabled')
p=Path(__file__).resolve().parent
for n in [18,19,20]:
 subprocess.run([sys.executable,str(p.parent/f'dim{n}/refinement_cyclic/verify.py')],check=True)
subprocess.run([sys.executable,str(p/'verify_algebra.py')],check=True)
print('PASS inherited full cyclic domains and saturated bounds in 18D-20D')
