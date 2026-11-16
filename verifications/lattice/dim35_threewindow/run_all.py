from pathlib import Path
import sys,subprocess,hashlib,json
from fractions import Fraction as F
if not __debug__:raise RuntimeError('Assertions required')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent
b=(p/'certificate.json').read_bytes()
assert hashlib.sha256(b).hexdigest()=='e6ba7b62cc3a42d5e4e28282b0e69d068536e34768a37b9d9f43c03e2ac367e8'
d=json.loads(b);assert d['n']==35 and F(d['safe_density'])==F('.000037664745') and 0<F(d['density'])<=F(d['safe_density'])
subprocess.run([sys.executable,str(p/'verify_independent.py'),str(p/'certificate.json'),'--bernstein-tail'],cwd=p,check=True)
print('COMPLETE 35D MULTI3 PASS <=0.000037664745',flush=True)
