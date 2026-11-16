from pathlib import Path
import subprocess,sys,json,hashlib
from fractions import Fraction as F
if not __debug__:raise RuntimeError('Assertions required; do not run with -O')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent
raw=(p/'certificate.json').read_bytes()
assert hashlib.sha256(raw).hexdigest()=='baddc0a1a28dc73f7e2cc9a779a2a2135616c1859aa595afa8bd73760efc14f4'
d=json.loads(raw);assert d['n']==33 and F(d.get('safe_density',d['density']))==F('0.000073827123') and 0<F(d['density'])<=F('0.000073827123')
subprocess.run([sys.executable,str(p/'verify_complete.py')],cwd=p,check=True)
print('COMPLETE 33D LATTICE CERTIFICATE PASS: density <=0.000073827123',flush=True)
