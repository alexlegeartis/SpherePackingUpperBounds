from pathlib import Path
import subprocess,sys,json,hashlib
from fractions import Fraction as F
if not __debug__:raise RuntimeError('Assertions required; do not run with -O')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent
raw=(p/'certificate.json').read_bytes()
assert hashlib.sha256(raw).hexdigest()=='766f55e618c395183f21c42b95b77cf3ba0211179904eaf4593d2ed98f07aa99'
d=json.loads(raw);assert d['n']==35 and F(d.get('safe_density',d['density']))==F('0.000038509364') and 0<F(d['density'])<=F('0.000038509364')
subprocess.run([sys.executable,str(p/'verify_complete.py')],cwd=p,check=True)
print('COMPLETE 35D LATTICE CERTIFICATE PASS: density <=0.000038509364',flush=True)
