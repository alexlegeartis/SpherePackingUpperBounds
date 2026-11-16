from pathlib import Path
import subprocess,sys,json,hashlib
from fractions import Fraction as F
if not __debug__:raise RuntimeError('Assertions required; do not run with -O')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent
raw=(p/'certificate.json').read_bytes()
assert hashlib.sha256(raw).hexdigest()=='3b9d68155a1d03c7e7e282012888845520559e05ed797818523d79f5f8873d78'
d=json.loads(raw);assert d['n']==20 and F(d.get('safe_density',d['density']))==F('0.005954630716') and 0<F(d['density'])<=F('0.005954630716')
subprocess.run([sys.executable,str(p/'verify_complete.py')],cwd=p,check=True)
print('COMPLETE 20D LATTICE CERTIFICATE PASS: density <=0.005954630716',flush=True)
