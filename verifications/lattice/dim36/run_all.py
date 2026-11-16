from pathlib import Path
import subprocess,sys,json,hashlib
from fractions import Fraction as F
if not __debug__:raise RuntimeError('Assertions required; do not run with -O')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent
raw=(p/'certificate.json').read_bytes()
assert hashlib.sha256(raw).hexdigest()=='66c5abd1ea57d5343cbc3a3f7d8dedd9699116afde305d219c06d5d5c3d6433c'
d=json.loads(raw);assert d['n']==36 and F(d.get('safe_density',d['density']))==F('0.000028235617') and 0<F(d['density'])<=F('0.000028235617')
subprocess.run([sys.executable,str(p/'verify_complete.py')],cwd=p,check=True)
print('COMPLETE 36D LATTICE CERTIFICATE PASS: density <=0.000028235617',flush=True)
