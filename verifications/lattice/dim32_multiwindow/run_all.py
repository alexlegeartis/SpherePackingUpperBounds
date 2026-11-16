from pathlib import Path
import sys,subprocess,hashlib,json
from fractions import Fraction as F
if not __debug__:raise RuntimeError('Assertions required')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent
b=(p/'certificate.json').read_bytes()
assert hashlib.sha256(b).hexdigest()=='7e79ec20cbcdda94e569193379b5197f230a24da0fff7e34c6b748ff3b1b7c87'
d=json.loads(b);assert d['n']==32 and F(d['safe_density'])==F('.000102834508') and 0<F(d['density'])<=F(d['safe_density'])
subprocess.run([sys.executable,str(p/'verify_independent.py'),str(p/'certificate.json'),'--bernstein-tail'],cwd=p,check=True)
print('COMPLETE 32D TWO-WINDOW PASS <=0.000102834508',flush=True)
