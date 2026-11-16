from pathlib import Path
import subprocess,sys,json,hashlib
from fractions import Fraction as F
if not __debug__:raise RuntimeError('Assertions required; do not run with -O')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent
raw=(p/'certificate.json').read_bytes()
assert hashlib.sha256(raw).hexdigest()=='dba4ba922b74fc1d1ba712a5161db80739ec8fce18af71963d4d7064b0eadb7c'
d=json.loads(raw);assert d['n']==27 and F(d.get('safe_density',d['density']))==F('0.000622938953') and 0<F(d['density'])<=F('0.000622938953')
subprocess.run([sys.executable,str(p/'verify_complete.py')],cwd=p,check=True)
print('COMPLETE 27D LATTICE CERTIFICATE PASS: density <=0.000622938953',flush=True)
