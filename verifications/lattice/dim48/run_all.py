from pathlib import Path
import sys,subprocess,hashlib,json,importlib.util
from fractions import Fraction as F
if not __debug__:raise RuntimeError('Assertions required; run without -O')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent
raw=(p/'certificate.json').read_bytes()
assert hashlib.sha256(raw).hexdigest()=='e8be45f13a6e3110846ae41c6bde16448b4389306478a68669932df6c236934e'
d=json.loads(raw)
assert d['n']==48 and F(d['gap'])==0 and F(d['safe_density'])==F('0.000000408175')
assert 0<F(d['density'])<=F(d['safe_density'])
flags=['--bernstein-tail']
if importlib.util.find_spec('flint') is not None:flags.append('--bernstein-fourier')
print('Full exact replay; Fourier method:', 'Bernstein (public python-flint)' if '--bernstein-fourier' in flags else 'SymPy root count',flush=True)
result=subprocess.run([sys.executable,str(p/'verify_independent.py'),str(p/'certificate.json'),*flags],cwd=p,check=True)
assert result.returncode==0
print('COMPLETE 48D PASS <=0.000000408175; verifier exit 0',flush=True)
