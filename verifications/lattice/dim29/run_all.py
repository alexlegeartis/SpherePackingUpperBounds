from pathlib import Path
import subprocess,sys,json,hashlib
from fractions import Fraction as F
if not __debug__:raise RuntimeError('Assertions required; do not run with -O')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent
raw=(p/'certificate.json').read_bytes()
assert hashlib.sha256(raw).hexdigest()=='a0ee675a01bc4ac521849b7d728044efeefbf10fc74ee59564e8c1a304864c57'
d=json.loads(raw);assert d['n']==29 and F(d.get('safe_density',d['density']))==F('0.000298025039') and 0<F(d['density'])<=F('0.000298025039')
subprocess.run([sys.executable,str(p/'verify_complete.py')],cwd=p,check=True)
print('COMPLETE 29D LATTICE CERTIFICATE PASS: density <=0.000298025039',flush=True)
