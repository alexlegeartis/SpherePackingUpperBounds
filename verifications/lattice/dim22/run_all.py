from pathlib import Path
import subprocess,sys,json,hashlib
from fractions import Fraction as F
if not __debug__:raise RuntimeError('Assertions required; do not run with -O')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent
raw=(p/'certificate.json').read_bytes()
assert hashlib.sha256(raw).hexdigest()=='6a4c69e9c17e1170c5521a49f9269ed216a5db5794eb261f6abd368fc71679bb'
d=json.loads(raw);assert d['n']==22 and F(d.get('safe_density',d['density']))==F('0.003482932334') and 0<F(d['density'])<=F('0.003482932334')
subprocess.run([sys.executable,str(p/'verify_complete.py')],cwd=p,check=True)
print('COMPLETE 22D LATTICE CERTIFICATE PASS: density <=0.003482932334',flush=True)
