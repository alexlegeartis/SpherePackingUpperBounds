"""Check every published count endpoint with the smaller half-open allowance."""
if not __debug__: raise RuntimeError('Verification must not run under python -O.')
from pathlib import Path
from fractions import Fraction as F
import json
from flint import fmpq as Q, fmpq_poly as Poly, arb, ctx
ctx.prec=256
root=Path(__file__).resolve().parent
s=json.loads((root/'witness.json').read_text());data=json.loads((root/'witness_psd_exact.json').read_text())
assert s['dimension']==5 and s['p']==4 and s['K']==12 and Q(s['R'])==Q(3,2)
M=[[[Q(v) for v in row] for row in block] for block in data['matrices']]
x=Poly([0,1]);z=4*x-5;ts=[Poly([1]),z]
for j in range(1,4):ts.append(2*z*ts[-1]-ts[-2])
q=sum((2*M[0][0][i+1]*ts[i] for i in range(5)),Poly([]))
for ell in range(13):
 o=int(ell==0)
 q+=sum((M[ell][i+o][j+o]*ts[i]*ts[j] for i in range(5) for j in range(5)),Poly([]))
ls=[Poly([1]),Q(5,2)-2*x]
for k in range(1,31):ls.append(((2*k+Q(5,2)-2*x)*ls[-1]-(k+Q(3,2))*ls[-2])/Q(k+1))
p=sum((Q(str(F(str(c))))*ls[k]/ls[k][0] for k,c in enumerate(s['x'])),Poly([]))
def ev(poly,t):
 v=arb(0)
 for c in reversed(poly.coeffs()):v=v*t+arb(c)
 return v
radii=list(map(F,s['count_radii']));weights=[F(str(w)) for w in s['count_weights']]
assert radii==list(map(F,['1.003012','1.008096','1.013258','1.017445']))
reports=[]
for radius in radii:
 r=Q(str(radius));X=arb.pi()*arb(r)**2
 allowance=sum((w for endpoint,w in zip(radii,weights) if radius<endpoint),F(0))
 value=ev(p,X)*(-X).exp()+arb(q(r))-arb(Q(str(allowance)))
 assert value<0,(str(radius),str(value))
 reports.append(dict(radius=str(radius),smaller_allowance=str(allowance),residual=str(value)))
(root/'count_endpoint_result.json').write_text(json.dumps(dict(verified=True,endpoints=reports),indent=2)+'\n')
print('Explicit half-open count endpoints: PASS',flush=True)
