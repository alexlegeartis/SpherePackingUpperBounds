"""Independent structure and infinite cyclic tail only. NOT a density proof."""
if not __debug__:
    raise RuntimeError('Assertions required')
from pathlib import Path
from fractions import Fraction as F
import json, hashlib
from flint import arb, fmpq, fmpq_poly, ctx
ctx.prec = 256
D = Path(__file__).resolve().parent
def rat(x):
    x = F(str(x))
    return fmpq(x.numerator, x.denominator)
def sturm(p):
    seq = [p, p.derivative()]
    while seq[-1].degree() > 0:
        q = -(seq[-2] % seq[-1])
        assert q
        seq.append(q)
    return seq
def variations(seq, x):
    values = [p(x) for p in seq]
    signs = [1 if z > 0 else -1 for z in values if z]
    return sum(a != b for a,b in zip(signs, signs[1:]))
source = json.loads((D/'angular.json').read_text())
fit = json.loads((D/'fit.json').read_text())
assert source['dimension'] == 9 if 'dimension' in source else True
c = list(map(rat, source['gegenbauer_coefficients']))
assert min(c) >= 0 and c[0] == 1
B = sum(c)
assert B == rat(source['B_exact'])
t = fmpq_poly([0,1])
G = [fmpq_poly([1]),t]
for k in range(1,len(c)-1):
    G.append(((2*k+7)*t*G[-1]-k*G[-2])/(k+7))
Q = sum((a*g for a,g in zip(c,G)),fmpq_poly())
assert list(Q.coeffs()) == list(map(rat,source['polynomial_coefficients']))
eps = fmpq(1,1000000)
q = Q + eps
assert q(-1) < 0 and q(fmpq(1,2)) < 0
sq = sturm(q)
assert variations(sq,-1) == variations(sq,fmpq(1,2))
sq = sturm(Q)
assert Q(fmpq(1,2)) < 0 and Q(1) > 0
assert variations(sq,fmpq(1,2))-variations(sq,1) == 1
p = list(map(rat,fit['parameters']))
nv=len(fit['v_knots']);nl=len(fit['log_knots']);assert len(p)==2+nv+nl
VK=list(map(rat,fit['v_knots']));LK=list(map(rat,fit['log_knots']))
assert all(1<k<=2 for k in VK+LK)
assert p[0]+p[1+nv]>=0
for k in [fmpq(1)]+VK:
 assert p[0]+sum((p[1+i] for i,x in enumerate(VK) if x<=k),fmpq(0))>=0
for k in sorted(set(VK+LK)):
 assert sum((k*p[1+i] for i,x in enumerate(VK) if x==k),fmpq(0))+sum((p[2+nv+i] for i,x in enumerate(LK) if x==k),fmpq(0))>=0
A=sum(p[:1+nv]);Bc=sum(p[1+nv:]);gamma=A+min(Bc,fmpq(0))/16
assert gamma>=7
phi16=arb(p[0])*15+arb(p[1+nv])*arb(16).log()
for i,k in enumerate(VK):phi16+=arb(p[1+i])*arb(16-k)
for i,k in enumerate(LK):phi16+=arb(p[2+nv+i])*(arb(16)/arb(k)).log()
assert phi16>0
tail1=arb(B)*(-phi16).exp()
tail2=arb(B)*((-3*arb(gamma)).exp()+(-phi16).exp())
assert tail1 < arb(eps) and tail2 < arb(eps)
out = dict(status='PASS', scope='Exact angular structure, all-scale ratio hypothesis, and u>=4 cyclic tail only. Compact cyclic, radial signs, pair and density rounding are NOT proved.',
    dimension=9, angular_coefficients_nonnegative=True, profile_g_monotone=True, one_sign_transition=True,
    cap_margin=str(eps), exact_final_linear_slope=str(gamma),
    tail_case1=str(tail1),tail_case2=str(tail2),
    hashes={f:hashlib.sha256((D/f).read_bytes()).hexdigest() for f in ['angular.json','fit.json']})
(D/'RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
