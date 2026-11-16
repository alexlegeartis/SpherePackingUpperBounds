from pathlib import Path
import sys,subprocess,hashlib,json
from fractions import Fraction as F
if not __debug__:raise RuntimeError('Assertions required')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent
b=(p/'certificate.json').read_bytes()
assert hashlib.sha256(b).hexdigest()=='4055c1f210d175c37e368c9652f881ff4bc62ba321f408a1909727e34ac4dc83'
d=json.loads(b);assert d['n']==32 and F(d['safe_density'])==F('.000101511777') and 0<F(d['density'])<=F(d['safe_density'])
subprocess.run([sys.executable,str(p/'verify_independent.py'),str(p/'certificate.json'),'--bernstein-tail'],cwd=p,check=True)
print('COMPLETE 32D REFINED PASS <=0.000101511777',flush=True)
