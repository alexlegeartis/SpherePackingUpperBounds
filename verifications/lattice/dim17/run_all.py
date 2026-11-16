from pathlib import Path
import subprocess,sys,json,hashlib
from fractions import Fraction
if not __debug__:raise RuntimeError('Assertions required; do not run with -O')
p=Path(__file__).resolve().parent
raw=(p/'certificate.json').read_bytes()
assert hashlib.sha256(raw).hexdigest()=='4eeb6ebef7e74c1f4a5d26af00f949411cfb57a8a468516cdf644778d346b626'
d=json.loads(raw);assert d['n']==17 and Fraction(d['density'])==Fraction('0.013832417384')
for script in ['verify_certificate.py','verify_independent.py']:
 subprocess.run([sys.executable,str(p/script),str(p/'certificate.json')],cwd=p,check=True)
print('COMPLETE 17D LATTICE CERTIFICATE PASS: density <=0.013832417384',flush=True)
