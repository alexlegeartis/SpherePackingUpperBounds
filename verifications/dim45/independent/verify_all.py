from pathlib import Path
import subprocess,sys
assert __debug__
p=Path(__file__).resolve().parent
for script in ['verify_angular_cyclic.py','verify.py']:
 subprocess.run([sys.executable,str(p/script)],check=True)
print('PASS complete independent45D new-angular and radial certificate')
