from pathlib import Path
import sys,subprocess,hashlib,json
from fractions import Fraction as F
if not __debug__:raise RuntimeError('Assertions required')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent
b=(p/'certificate.json').read_bytes()
assert hashlib.sha256(b).hexdigest()=='ee54744cac8fdb49c1764ebf2bde57e0b01a0cb967992282ab64779d0cc1ee1a'
d=json.loads(b);assert d['n']==37 and F(d['safe_density'])==F('.000020580470') and 0<F(d['density'])<=F(d['safe_density'])
subprocess.run([sys.executable,str(p/'verify_independent.py'),str(p/'certificate.json'),'--bernstein-tail'],cwd=p,check=True)
print('COMPLETE 37D PASS <=0.000020580470',flush=True)
