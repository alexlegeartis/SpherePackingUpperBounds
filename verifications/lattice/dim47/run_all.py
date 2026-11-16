from pathlib import Path
import sys,subprocess,hashlib,json,importlib.util
from fractions import Fraction as F
if not __debug__:raise RuntimeError('Assertions required; run without -O')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent
raw=(p/'certificate.json').read_bytes()
assert hashlib.sha256(raw).hexdigest()=='d1ea47251eb7d77ab1a7d3b119be1e8ae6515a0cd25db8c511e41b1495031ef8'
d=json.loads(raw)
assert d['n']==47 and F(d['gap'])==0 and F(d['safe_density'])==F('0.000000591010')
assert 0<F(d['density'])<=F(d['safe_density'])
flags=['--bernstein-tail']
if importlib.util.find_spec('flint') is not None:flags.append('--bernstein-fourier')
print('Full exact replay; Fourier method:', 'Bernstein (public python-flint)' if '--bernstein-fourier' in flags else 'SymPy root count',flush=True)
result=subprocess.run([sys.executable,str(p/'verify_independent.py'),str(p/'certificate.json'),*flags],cwd=p,check=True)
assert result.returncode==0
print('COMPLETE 47D PASS <=0.000000591010; verifier exit 0',flush=True)
