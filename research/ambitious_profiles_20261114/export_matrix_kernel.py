"""Exact PSD reconstruction, origin saturation and cap separation for matrix kernels."""
import argparse,json,sys,math
from pathlib import Path
from fractions import Fraction as F
from flint import fmpq as Q,fmpq_poly as P,arb,arb_poly,ctx
from common import ROOT,normalized_gegenbauer,rational
if not __debug__:raise RuntimeError('Assertions required')
sys.set_int_max_str_digits(0);ctx.prec=512
ap=argparse.ArgumentParser();ap.add_argument('candidate');ap.add_argument('--epsilon',default='1/10000')
args=ap.parse_args();file=Path(args.candidate).resolve();D=file.parent
r=json.loads(file.read_text());n=r['dimension'];eps=rational(args.epsilon)
assert r.get('fixed_slow_angular'), 'This exporter freezes the slow angular polynomial'
ref=json.loads((ROOT/r['reference_angular']).read_text())
bell=list(map(Q,ref['gegenbauer_coefficients']));B=Q(ref['B_exact']);degree=len(bell)-1
active=[k for k,z in enumerate(bell) if z]
qfloat=lambda z:Q(str(F(float(z))))
raw=[[[qfloat(z) for z in row] for row in matrix] for matrix in r['blocks']]
for k in active:raw[k][0][0]=qfloat(r['origin_block'][1][1])
for matrix in raw:
    assert matrix[0][1]==matrix[1][0]
delta=3*eps
while True:
    mats=[]
    for k,m in enumerate(raw):
        a,c,d=m[0][0],m[0][1],m[1][1]
        if k in active:a+=delta;d+=delta
        if k==0:a-=eps;c-=eps;d-=eps
        mats.append([[a,c],[c,d]])
    if all(mats[k][0][0]>0 and mats[k][1][1]>0 and mats[k][0][0]*mats[k][1][1]-mats[k][0][1]**2>0 for k in active):break
    delta*=2
    assert delta<1
cross=[qfloat(r['origin_block'][0][1])-delta,qfloat(r['origin_block'][0][2])-delta]
a,c,d=mats[0][0][0],mats[0][0][1],mats[0][1][1]
det=a*d-c*c
origin=(d*cross[0]**2-2*c*cross[0]*cross[1]+a*cross[1]**2)/det
assert origin>0
mats=[[[z/origin for z in row] for row in matrix] for matrix in mats]
cross=[z/origin for z in cross]
assert cross[0]<0
G=normalized_gegenbauer(n,degree)
polys=[]
for i,j in [(0,0),(0,1),(1,1)]:
    polys.append(sum((bell[k]*mats[k][i][j]*G[k] for k in active),P([])))
refpoly=P(list(map(Q,ref['polynomial_coefficients'])))
scale=(qfloat(r['origin_block'][1][1])+delta)/origin
subtract=eps/origin
assert polys[0]==scale*refpoly-subtract
strict=subtract/10
bad=[]
for name,poly in [('slow_cap',polys[0]),('mixed_cap',polys[0]+polys[1])]:
    apoly=arb_poly([arb(c) for c in poly.coeffs()])
    candidates=[(-1.,apoly(arb(-1))),(.5001,apoly(arb(Q(5001,10000))))]
    for z,m in poly.derivative().numer().complex_roots():
        if z.imag.contains(0) and z.real>arb(-1) and z.real<arb(Q(5001,10000)):
            candidates.append((float(z.real.mid()),apoly(z.real)))
    best=max(candidates,key=lambda z:float(z[1].mid()))
    print(name,'maximum',str(best[1]),'point',best[0],flush=True)
    if not best[1]<-arb(strict):
        bad += [x for x,y in candidates if not y<-arb(strict)]
if bad:
    cutsfile=D/'cuts.json';cuts=json.loads(cutsfile.read_text()) if cutsfile.exists() else {}
    cuts['cap']=sorted(set(cuts.get('cap',[])+bad))
    cutsfile.write_text(json.dumps(cuts,indent=2)+'\n')
    (D/'kernel_cap_rejection.json').write_text(json.dumps({'cap_points':bad},indent=2)+'\n')
    raise RuntimeError('Exact cap separation requires another SDP solve')
# Sturm excludes hidden positive cap intervals; root maxima above only choose cuts.
def variation(poly,point):
    seq=[poly,poly.derivative()]
    while seq[-1].degree()>0:
        rem=-(seq[-2]%seq[-1]);assert rem
        num=rem.numer();seq.append(P(num/abs(num.content())))
    vals=[p(point) for p in seq];sg=[1 if x>0 else -1 for x in vals if x]
    return sum(x!=y for x,y in zip(sg,sg[1:]))
for poly in [polys[0]+strict,polys[0]+polys[1]+strict]:
    assert poly(-1)<0 and poly(Q(5001,10000))<0
    assert variation(poly,Q(-1))==variation(poly,Q(5001,10000))

Babs=[sum((bell[k]*abs(mats[k][i][j]) for k in active),Q(0)) for i,j in [(0,0),(0,1),(1,1)]]
slow,fast=rational(r['slow_exact']),rational(r['fast_exact']);gamma=fast-slow
assert slow>=40 and gamma>0
eta=subtract
Tmax=Q(2)
while not 6*arb(Babs[1])*(-arb(gamma)*(arb(Tmax)-1)).exp()+3*arb(Babs[2])*(-2*arb(gamma)*(arb(Tmax)-1)).exp()<arb(eta)/2:
    Tmax+=1
U=Q(5);cut=Q(5001,10000)
tail=(2*arb(Babs[1]+Babs[2])*(-arb(gamma)*(U*U-1)).exp()
      +arb(Babs[0]+Babs[1])*(-arb(slow)*(2*cut*U-1)).exp()
      +arb(Babs[0]+2*Babs[1]+Babs[2])*(-arb(slow)*(U*U-1)).exp())
assert tail<arb(strict)/2
out={'dimension':n,'angular_degree':degree,'s_exact':r['s_exact'],
     'slow_exact':str(slow),'fast_exact':str(fast),'reference_angular':r['reference_angular'],
     'reference_B_exact':str(B),'reference_coefficients':list(map(str,bell)),
     'blocks':[[[str(z) for z in row] for row in matrix] for matrix in mats],
     'cross_vector_exact':list(map(str,cross)),'origin_exact':'1',
     'origin_before_normalization_exact':str(origin),
     'origin_schur_saturated_exact':True,'PSD_regularizer_exact':str(delta),'constant_subtraction_exact':str(eps),
     'slow_polynomial_reference_scale_exact':str(scale),'slow_polynomial_constant_subtraction_exact':str(subtract),
     'angular_polynomials':[list(map(str,p.coeffs())) for p in polys],
     'angular_absolute_bounds':list(map(str,Babs)),
     'cap_margin_exact':str(strict),'large_T_margin_exact':str(eta),
     'compact_T_upper_exact':str(Tmax),'compact_u_upper_exact':str(U),
     'cap_cutoff_exact':str(cut),'certified':False,
     'scope':'Exact PSD/origin/cap and analytic tails; full compact triangle proof remains required'}
(D/'kernel.json').write_text(json.dumps(out,indent=2)+'\n')
print('EXACT KERNEL PSD/CAP/TAIL PASS','origin',float(origin),'Tmax',str(Tmax),'epsilon',str(eps),flush=True)
