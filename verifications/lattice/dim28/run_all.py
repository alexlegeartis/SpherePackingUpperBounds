from pathlib import Path
import subprocess,sys,json,hashlib
from fractions import Fraction as F
if not __debug__:raise RuntimeError('Assertions required; do not run with -O')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent
raw=(p/'certificate.json').read_bytes()
assert hashlib.sha256(raw).hexdigest()=='102f608dd7d17a26d583e9d6918a23520077369d6a54db747dc10430884b3418'
d=json.loads(raw);assert d['n']==28 and F(d.get('safe_density',d['density']))==F('0.000429639101') and 0<F(d['density'])<=F('0.000429639101')
subprocess.run([sys.executable,str(p/'verify_complete.py')],cwd=p,check=True)
print('COMPLETE 28D LATTICE CERTIFICATE PASS: density <=0.000429639101',flush=True)
