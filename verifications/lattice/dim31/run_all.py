from pathlib import Path
import subprocess,sys,json,hashlib
from fractions import Fraction as F
if not __debug__:raise RuntimeError('Assertions required; do not run with -O')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent
raw=(p/'certificate.json').read_bytes()
assert hashlib.sha256(raw).hexdigest()=='b29a49db41b954edee4fcc93af836f87a9d537e998b3a60fde91b01399119a49'
d=json.loads(raw);assert d['n']==31 and F(d.get('safe_density',d['density']))==F('0.000146087217') and 0<F(d['density'])<=F('0.000146087217')
subprocess.run([sys.executable,str(p/'verify_complete.py')],cwd=p,check=True)
print('COMPLETE 31D LATTICE CERTIFICATE PASS: density <=0.000146087217',flush=True)
