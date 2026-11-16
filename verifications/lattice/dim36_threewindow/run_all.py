from pathlib import Path
import sys,subprocess,hashlib,json
from fractions import Fraction as F
if not __debug__:raise RuntimeError('Assertions required')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent
b=(p/'certificate.json').read_bytes()
assert hashlib.sha256(b).hexdigest()=='08534334526e55c9d38b99fd9b1643ff8f0104c9c0b1269375c069e8ff8d60a3'
d=json.loads(b);assert d['n']==36 and F(d['safe_density'])==F('.000027428605') and 0<F(d['density'])<=F(d['safe_density'])
subprocess.run([sys.executable,str(p/'verify_complete.py')],cwd=p,check=True)
print('COMPLETE 36D MULTI3 PASS <=0.000027428605',flush=True)
