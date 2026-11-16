from pathlib import Path
import subprocess,sys,json,hashlib
from fractions import Fraction as F
if not __debug__:raise RuntimeError('Assertions required; do not run with -O')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent
raw=(p/'certificate.json').read_bytes()
assert hashlib.sha256(raw).hexdigest()=='014502c998b3e16dccf1c7c298b0f250b3df3c24f72f3087c5acf1e95c38e065'
d=json.loads(raw);assert d['n']==32 and F(d['safe_density'])==F('.000104134823') and 0<F(d['density'])<=F('.000104134823')
subprocess.run([sys.executable,str(p/'verify_complete.py')],cwd=p,check=True)
print('COMPLETE 32D LATTICE CERTIFICATE PASS: density <=.000104134823',flush=True)
