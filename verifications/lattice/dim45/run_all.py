from pathlib import Path
import sys,subprocess,hashlib,json,importlib.util
from fractions import Fraction as F
if not __debug__:raise RuntimeError('Assertions required; run without -O')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent
raw=(p/'certificate.json').read_bytes()
assert hashlib.sha256(raw).hexdigest()=='e32c22fea52e563189a2aebb6ba637ddb8c5c45ce862cbd5ba8b1b1facd9f49d'
d=json.loads(raw)
assert d['n']==45 and F(d['gap'])==0 and F(d['safe_density'])==F('0.000001263630')
assert 0<F(d['density'])<=F(d['safe_density'])
flags=['--bernstein-tail']
if importlib.util.find_spec('flint') is not None:flags.append('--bernstein-fourier')
print('Full exact replay; Fourier method:', 'Bernstein (public python-flint)' if '--bernstein-fourier' in flags else 'SymPy root count',flush=True)
result=subprocess.run([sys.executable,str(p/'verify_independent.py'),str(p/'certificate.json'),*flags],cwd=p,check=True)
assert result.returncode==0
print('COMPLETE 45D PASS <=0.000001263630; verifier exit 0',flush=True)
