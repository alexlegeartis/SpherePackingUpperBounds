from pathlib import Path
import subprocess,sys
assert __debug__
p=Path(__file__).resolve().parent
for s in ['algebra_pair.py','domains.py']:
 subprocess.run([sys.executable,str(p/s)],check=True)
print('PASS independent full cyclic certificate: density <0.018250915, Schur tightened <0.018204096')
