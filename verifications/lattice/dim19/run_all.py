from pathlib import Path
import subprocess,sys,json,hashlib
from fractions import Fraction
if not __debug__:raise RuntimeError('Assertions required; do not run with -O')
p=Path(__file__).resolve().parent
raw=(p/'certificate.json').read_bytes()
assert hashlib.sha256(raw).hexdigest()=='7196d4c4c1b0e33896413ca2ef064cef6c0ffecd77049b49532010964ec07447'
d=json.loads(raw);assert d['n']==19 and Fraction(d['density'])==Fraction('0.007808430327')
for script in ['verify_certificate.py','verify_independent.py']:
 subprocess.run([sys.executable,str(p/script),str(p/'certificate.json')],cwd=p,check=True)
print('COMPLETE 19D LATTICE CERTIFICATE PASS: density <=0.007808430327',flush=True)
