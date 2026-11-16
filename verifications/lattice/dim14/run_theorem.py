if not __debug__:
    raise RuntimeError('Verification requires Python without -O')
from pathlib import Path
import subprocess,sys
P=Path(__file__).resolve().parent
import json
from fractions import Fraction as F
prior=json.loads((P/'prior_certificate.json').read_text())
assert prior['n']==14 and F(prior['claims']['density'])<=F('.03525')
subprocess.run([sys.executable,str(P.parent/'prior'/'taylor.py'),str(P/'prior_certificate.json')],cwd=P,check=True)
subprocess.run([sys.executable,str(P/'verify_theorem.py')],cwd=P,check=True)
