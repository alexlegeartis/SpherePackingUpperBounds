from pathlib import Path
import subprocess,sys,json,hashlib
from fractions import Fraction as F
if not __debug__:raise RuntimeError('Assertions required; do not run with -O')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent
raw=(p/'certificate.json').read_bytes()
assert hashlib.sha256(raw).hexdigest()=='f83d66a3746811c8887a0095a3df112ff0b322be24981634c02e493ae022778a'
d=json.loads(raw);assert d['n']==26 and F(d['safe_density'])==F('.000907209997') and 0<F(d['density'])<=F('.000907209997')
subprocess.run([sys.executable,str(p/'verify_complete.py')],cwd=p,check=True)
print('COMPLETE 26D LATTICE CERTIFICATE PASS: density <=.000907209997',flush=True)
