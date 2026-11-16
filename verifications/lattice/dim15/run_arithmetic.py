"""Portable exact15D prior, determinant, and arithmetic replay. Python3+SymPy."""
import subprocess,sys
from pathlib import Path
if not __debug__:raise RuntimeError('Assertions required; do not use -O')
here=Path(__file__).resolve().parent
for argv in [[str(here.parent/'prior'/'taylor.py'),str(here/'prior_certificate.json')],[str(here/'verify_determinant.py')],[str(here/'verify_arithmetic.py')]]:
 subprocess.run([sys.executable,*argv],cwd=here,check=True)
print('PASS:15D prior and arithmetic theorem dependencies')
