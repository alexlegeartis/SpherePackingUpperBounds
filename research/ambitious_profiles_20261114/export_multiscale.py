"""Rational PSD/origin reconstruction and Bernstein long-edge cap proof."""
import argparse,json,sys,math
from pathlib import Path
from fractions import Fraction as F
from flint import fmpq as Q,fmpq_poly as P,fmpq_mat,arb,arb_poly,ctx
from common import ROOT,normalized_gegenbauer,rational
if not __debug__:raise RuntimeError('Assertions required')
sys.set_int_max_str_digits(0);ctx.prec=512
ap=argparse.ArgumentParser();ap.add_argument('candidate');ap.add_argument('--epsilon',default='1/100000')
args=ap.parse_args();file=Path(args.candidate).resolve();D=file.parent;r=json.loads(file.read_text())
n=r['dimension'];rates=list(map(Q,r['rates_exact']));m=len(rates);powers=r['powers'];gamma=rates[1]-rates[0]
assert rates[0]>=40 and all(gamma*k==b-rates[0] for k,b in zip(powers,rates))
ref=json.loads((ROOT/r['reference_angular']).read_text());bell=list(map(Q,ref['gegenbauer_coefficients']));B=Q(ref['B_exact'])
d=len(bell)-1;active=[k for k,z in enumerate(bell) if z];indices=[(i,j) for i in range(m) for j in range(i,m)]
f=lambda z:Q(str(F(float(z))))
raw=[[[f(z) for z in row] for row in matrix] for matrix in r['blocks']]
common=f(r['origin_block'][1][1]);eps=rational(args.epsilon);delta=(m+1)*eps
for k in active:raw[k][0][0]=common
def positive_definite(A):
    W=[row[:] for row in A]
    for j in range(len(A)):
        pivot=W[j][j]
        if pivot<=0:return False
        for i in range(j+1,len(A)):
            for k in range(i,len(A)):
                W[k][i]=W[i][k]=W[i][k]-W[i][j]*W[k][j]/pivot
    return True
while True:
    matrices=[]
    for ell,A in enumerate(raw):
        M=[[A[i][j]+(delta if i==j and ell in active else 0)-(eps if ell==0 else 0) for j in range(m)] for i in range(m)]
        matrices.append(M)
    if all(positive_definite(matrices[k]) for k in active):break
    delta*=2;assert delta<1
cross=[f(r['origin_block'][0][i+1])-delta for i in range(m)]
inverse=fmpq_mat(matrices[0]).inv()
origin=sum((cross[i]*inverse[i,j]*cross[j] for i in range(m) for j in range(m)),Q(0))
assert origin>0
matrices=[[[z/origin for z in row] for row in A] for A in matrices];cross=[z/origin for z in cross]
G=normalized_gegenbauer(n,d)
polys=[sum((bell[ell]*matrices[ell][i][j]*G[ell] for ell in active),P([])) for i,j in indices]
refpoly=P(list(map(Q,ref['polynomial_coefficients'])))
scale=(common+delta)/origin;eta=eps/origin;strict=eta/10
assert polys[0]==scale*refpoly-eta
K=max(powers);cap_polys=[]
for j in range(K+1):
    poly=polys[indices.index((0,0))]
    for i in range(1,m):
        if j>=powers[i]:poly+=Q(math.comb(j,powers[i]),math.comb(K,powers[i]))*polys[indices.index((0,i))]
    cap_polys.append(poly)
bad=[];cut=Q(5001,10000)
for j,poly in enumerate(cap_polys):
    AP=arb_poly([arb(z) for z in poly.coeffs()]);points=[(-1.,AP(arb(-1))),(.5001,AP(arb(cut)))]
    for z,mul in poly.derivative().numer().complex_roots():
        if z.imag.contains(0) and z.real>arb(-1) and z.real<arb(cut):points.append((float(z.real.mid()),AP(z.real)))
    best=max(points,key=lambda z:float(z[1].mid()));print('CAP',j,'maximum',str(best[1]),'at',best[0],flush=True)
    bad.extend(x for x,y in points if not y<-arb(strict))
if bad:
    cp=D/'cuts.json';cuts=json.loads(cp.read_text()) if cp.exists() else {}
    cuts['cap']=sorted(set(cuts.get('cap',[])+bad));cp.write_text(json.dumps(cuts,indent=2)+'\n')
    (D/'cap_rejection.json').write_text(json.dumps({'counterexample_candidates':bad},indent=2)+'\n')
    raise RuntimeError('Cap separation requires another solve')
def variations(poly,x):
    seq=[poly,poly.derivative()]
    while seq[-1].degree()>0:
        rem=-(seq[-2]%seq[-1]);assert rem
        z=rem.numer();seq.append(P(z/abs(z.content())))
    vals=[p(x) for p in seq];sg=[1 if z>0 else -1 for z in vals if z]
    return sum(x!=y for x,y in zip(sg,sg[1:]))
for poly in cap_polys:
    poly+=strict;assert poly(-1)<0 and poly(cut)<0 and variations(poly,Q(-1))==variations(poly,cut)
Babs=[sum((bell[k]*abs(matrices[k][i][j]) for k in active),Q(0)) for i,j in indices]
def bij(i,j):return Babs[indices.index((min(i,j),max(i,j))) ]
slowfast=sum((bij(0,i) for i in range(1,m)),Q(0))
fastfast=sum((bij(i,j) for i in range(1,m) for j in range(1,m)),Q(0))
colfast=sum((bij(i,j) for i in range(m) for j in range(1,m)),Q(0))
col0=sum((bij(i,0) for i in range(m)),Q(0));total=sum((bij(i,j) for i in range(m) for j in range(m)),Q(0))
Tmax=Q(2)
while not 6*arb(slowfast)*(-arb(gamma)*(Tmax-1)).exp()+3*arb(fastfast)*(-2*arb(gamma)*(Tmax-1)).exp()<arb(eta)/2:Tmax+=1
Umax=Q(5)
tail=2*arb(colfast)*(-arb(gamma)*(Umax*Umax-1)).exp()+arb(col0)*(-arb(rates[0])*(2*cut*Umax-1)).exp()+arb(total)*(-arb(rates[0])*(Umax*Umax-1)).exp()
assert tail<arb(strict)/2 and cross[0]<0
out={'dimension':n,'features':m,'angular_degree':d,'s_exact':r['s_exact'],
     'rates_exact':list(map(str,rates)),'powers':powers,'reference_angular':r['reference_angular'],
     'reference_B_exact':str(B),'reference_coefficients':list(map(str,bell)),
     'blocks':[[[str(z) for z in row] for row in A] for A in matrices],
     'cross_vector_exact':list(map(str,cross)),'origin_exact':'1','origin_before_normalization_exact':str(origin),
     'PSD_regularizer_exact':str(delta),'constant_subtraction_exact':str(eps),
     'slow_polynomial_reference_scale_exact':str(scale),'slow_polynomial_constant_subtraction_exact':str(eta),
     'angular_indices':indices,'angular_polynomials':[list(map(str,p.coeffs())) for p in polys],
     'cap_Bernstein_polynomials':[list(map(str,p.coeffs())) for p in cap_polys],
     'angular_absolute_bounds':list(map(str,Babs)),'cap_margin_exact':str(strict),'large_T_margin_exact':str(eta),
     'compact_T_upper_exact':str(Tmax),'compact_u_upper_exact':str(Umax),'cap_cutoff_exact':str(cut),
     'certified':False,'scope':'PSD, saturated origin, Bernstein caps and analytic tails only; compact triangle proof required'}
(D/'kernel.json').write_text(json.dumps(out,indent=2)+'\n')
print('EXACT MULTISCALE PSD/CAP/TAIL PASS','origin',float(origin),'Tmax',str(Tmax),flush=True)
