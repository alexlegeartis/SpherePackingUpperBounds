from pathlib import Path
import subprocess,sys,json,hashlib
from fractions import Fraction as F
if not __debug__:raise RuntimeError('Assertions required; do not run with -O')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent
raw=(p/'certificate.json').read_bytes()
assert hashlib.sha256(raw).hexdigest()=='12f13a429a21813a0a28f0f560b23e27d5e603727c977d44e666d51f4438e369'
d=json.loads(raw);assert d['n']==25 and F(d.get('safe_density',d['density']))==F('0.001325232433') and 0<F(d['density'])<=F('0.001325232433')
subprocess.run([sys.executable,str(p/'verify_complete.py')],cwd=p,check=True)
print('COMPLETE 25D LATTICE CERTIFICATE PASS: density <=0.001325232433',flush=True)
