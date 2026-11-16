from pathlib import Path
import subprocess,sys,json,hashlib
from fractions import Fraction as F
if not __debug__:raise RuntimeError('Assertions required; do not run with -O')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent
raw=(p/'certificate.json').read_bytes()
assert hashlib.sha256(raw).hexdigest()=='a82331fdb260af9bc90de5ab0df306c718f8123b71cab81f620c80122d787b9f'
d=json.loads(raw);assert d['n']==34 and F(d.get('safe_density',d['density']))==F('0.000053148401') and 0<F(d['density'])<=F('0.000053148401')
subprocess.run([sys.executable,str(p/'verify_complete.py')],cwd=p,check=True)
print('COMPLETE 34D LATTICE CERTIFICATE PASS: density <=0.000053148401',flush=True)
