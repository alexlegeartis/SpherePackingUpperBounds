from pathlib import Path
import subprocess,sys,json,hashlib
from fractions import Fraction as F
if not __debug__:raise RuntimeError('Assertions required; do not run with -O')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent
raw=(p/'certificate.json').read_bytes()
assert hashlib.sha256(raw).hexdigest()=='ed71fcdae1b7a4507449d951a8fd1a0eacbb7a4440d7e2e883ce5bf32faea823'
d=json.loads(raw);assert d['n']==21 and F(d.get('safe_density',d['density']))==F('0.004557959175') and 0<F(d['density'])<=F('0.004557959175')
subprocess.run([sys.executable,str(p/'verify_complete.py')],cwd=p,check=True)
print('COMPLETE 21D LATTICE CERTIFICATE PASS: density <=0.004557959175',flush=True)
