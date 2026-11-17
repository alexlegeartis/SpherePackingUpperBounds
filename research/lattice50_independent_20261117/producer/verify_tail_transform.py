from pathlib import Path
import json,math,sys
from flint import fmpq as Q,fmpq_poly as P,ctx
if not __debug__:raise RuntimeError('Assertions must remain enabled')
sys.set_int_max_str_digits(0);ctx.prec=1024;D=Path(__file__).parent;r=json.load(open(D/'repaired.json'));p=P(list(map(Q,r['P'])));q=P(list(map(Q,r['Q'])));deg=p.degree()
# Independent closed-form Gaussian monomial transform.
reconstructed=[]
for j in range(deg+1):
 reconstructed.append(sum((p[k]*math.factorial(k)*(-1)**j*math.prod(range(25+j,25+k))/(math.factorial(k-j)*math.factorial(j))for k in range(j,deg+1)),Q(0)))
assert P(reconstructed)==q and q(0)==1
roots=p.complex_roots();assert sum(m for z,m in roots)==deg;assert all(not z.imag.contains(0)or z.real<180 for z,m in roots);assert p(Q(180))<0
out=dict(status='INDEPENDENT EXACT TRANSFORM AND GLOBAL SPATIAL TAIL PASS',precision=1024,cut_x='180',degree=deg,total_root_multiplicity=sum(m for z,m in roots),possibly_real_root_count=sum(m for z,m in roots if z.imag.contains(0)),all_possibly_real_roots_below_cut=True,negative_at_cut=True)
(D/'independent_tail_transform_result.json').write_text(json.dumps(out,indent=2));print(out,flush=True)
