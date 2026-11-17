"""Independent closed-form harmonic reconstruction and monotone ratio enclosures."""
import argparse,json,math,sys,time
from pathlib import Path
from fractions import Fraction
from flint import fmpq as Q,fmpq_poly as P,arb,arb_poly,ctx
if not __debug__:raise RuntimeError('Assertions must remain enabled')
sys.set_int_max_str_digits(0);ctx.prec=512
ap=argparse.ArgumentParser();ap.add_argument('package');args=ap.parse_args();D=Path(args.package).resolve()
r=json.loads((D/'angular.json').read_text());cert=json.loads((D/'certificate.json').read_text())
qf=lambda z:Q(str(Fraction(str(z))))
n,d=r['dimension'],r['degree'];assert n>=3 and n==cert['dimension']
a,b=qf(r['profile_a_exact']),qf(r['profile_b_exact'])
assert a+b>0 and b>0 and a==qf(cert['a_exact']) and b==qf(cert['b_exact'])
poly=P(list(map(Q,r['polynomial_coefficients'])));cs=list(map(Q,r['gegenbauer_coefficients']))
delta,cut=qf(r['negative_cap_delta']),qf(r['cutoff']);B=qf(r['B_exact'])
assert len(cs)==d+1 and min(cs)>=0 and cs[0]==1 and Q(1,2)<cut<1
assert delta>0 and B>0 and d==poly.degree()
def rise(x,k):return math.prod((x+j for j in range(k)),start=Q(1))
lam=Q(n-2,2);rebuilt=P([])
for ell,coefficient in enumerate(cs):
    g=[Q(0)]*(ell+1);normal=rise(2*lam,ell)/math.factorial(ell)
    for j in range(ell//2+1):
        g[ell-2*j]=(-1)**j*rise(lam,ell-j)*2**(ell-2*j)/(math.factorial(j)*math.factorial(ell-2*j)*normal)
    gp=P(g);assert gp(1)==1;rebuilt+=coefficient*gp
assert rebuilt==poly and sum(cs)==B==poly(1)
moment=Q(0)
for j in range(poly.degree()//2+1):moment+=poly[2*j]*rise(Q(1,2),j)/rise(Q(n,2),j)
assert moment==1
def variations(p,x):
    seq=[p,p.derivative()]
    while seq[-1].degree()>0:
        rem=-(seq[-2]%seq[-1]);assert rem
        z=rem.numer();seq.append(P(z/abs(z.content())))
    vals=[z(x) for z in seq];sg=[1 if z>0 else -1 for z in vals if z]
    return sum(x!=y for x,y in zip(sg,sg[1:]))
assert poly(-1)<0 and poly(cut)<0 and poly(1)>0
assert variations(poly,Q(-1))==variations(poly,cut)
assert variations(poly,cut)-variations(poly,Q(1))==1
# Strict negative cap by exact Bernstein coefficients, separate from the Sturm proof.
cap=poly+delta;stack=[(Q(-1),Q(1,2),0)];cap_leaves=cap_depth=0
while stack:
    lo,hi,depth=stack.pop();translated=cap(P([lo,hi-lo]))
    bs=[sum((translated[j]*Q(math.comb(i,j),math.comb(d,j)) for j in range(i+1)),Q(0)) for i in range(d+1)]
    if max(bs)<0:cap_leaves+=1;cap_depth=max(cap_depth,depth);continue
    assert depth<35
    mid=(lo+hi)/2;stack.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
U=Q(r['tail_u']);power=max(Q(0),-a)
while True:
    assert U*U>power/b
    tail=arb(B)*(arb(4)**arb(power)*(-arb(b)*(2*cut*U-1)).exp()+arb(U)**arb(2*power)*(-arb(b)*(U*U-1)).exp())
    if tail<arb(delta)/2:break
    U*=2

assert U==8
res=dict(status='EXACT ANGULAR STRUCTURE AND INFINITE TAIL PASS',dimension=n,precision=ctx.prec,cap_leaves=cap_leaves,cap_depth=cap_depth,negative_cap=str(delta),tail_u=str(U),tail_bound=str(tail),B=str(B))
(D/'angular_structure_result.json').write_text(json.dumps(res,indent=2));print(res,flush=True)
