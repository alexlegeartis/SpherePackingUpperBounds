"""Exact alternative to real-root isolation for the infinite physical tail."""
import sys;sys.set_int_max_str_digits(0)
from pathlib import Path
import json,hashlib
from flint import fmpq_poly as Poly,fmpq as Q
if not __debug__:raise RuntimeError('Do not use python -O')
D=Path(__file__).resolve().parent;source=D/'input/root_repair/repaired.json';r=json.loads(source.read_text());coeff=list(map(Q,r['P']));p=Poly(coeff);shifted=Poly([0]);x=Poly([180,1])
for a in coeff[::-1]:shifted=shifted*x+a
assert shifted.degree()==p.degree() and all(c<0 for c in shifted)
# P(180+t) is a polynomial with every coefficient strictly negative; hence P(x)<0 for x>=180.
out=dict(status='EXACT TRANSLATED-COEFFICIENT TAIL PASS',n=r['dimension'],negative_from=180,degree=p.degree(),all_shifted_coefficients_strictly_negative=True,input_sha256=hashlib.sha256(source.read_bytes()).hexdigest());(D/'exact_tail_result.json').write_text(json.dumps(out,indent=2));print(out)
