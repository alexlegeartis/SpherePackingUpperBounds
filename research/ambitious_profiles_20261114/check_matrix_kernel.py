"""Exact matrix reconstruction and complete compact three-variable cyclic proof.
The unchanged reference Gaussian angular certificate must also be replayed.
"""
import argparse,json,math,sys,time
from pathlib import Path
from fractions import Fraction as F
from flint import fmpq as Q,fmpq_poly as P,arb,arb_poly,ctx
from common import ROOT,normalized_gegenbauer,rational
if not __debug__:raise RuntimeError('Assertions required')
sys.set_int_max_str_digits(0);ctx.prec=512
ap=argparse.ArgumentParser();ap.add_argument('kernel');ap.add_argument('--budget',type=int,default=500000)
ap.add_argument('--maxdepth',type=int,default=32);args=ap.parse_args()
file=Path(args.kernel).resolve();D=file.parent;r=json.loads(file.read_text())
n,d=r['dimension'],r['angular_degree'];bell=list(map(Q,r['reference_coefficients']));B=Q(r['reference_B_exact'])
mats=[[[Q(z) for z in row] for row in m] for m in r['blocks']]
cross=list(map(Q,r['cross_vector_exact']));slow,fast=Q(r['slow_exact']),Q(r['fast_exact']);gamma=fast-slow
ref=json.loads((ROOT/r['reference_angular']).read_text());refpoly=P(list(map(Q,ref['polynomial_coefficients'])))
assert bell==list(map(Q,ref['gegenbauer_coefficients'])) and sum(bell)==B
G=normalized_gegenbauer(n,d)
for k,m in enumerate(mats):
    a,c,e=m[0][0],m[0][1],m[1][1];assert c==m[1][0]
    if not bell[k]:assert a==c==e==0
    else:assert a>0 and e>0 and a*e-c*c>0
a,c,e=mats[0][0][0],mats[0][0][1],mats[0][1][1]
assert (e*cross[0]**2-2*c*cross[0]*cross[1]+a*cross[1]**2)/(a*e-c*c)==1
polys=[P(list(map(Q,p))) for p in r['angular_polynomials']]
for poly,(i,j) in zip(polys,[(0,0),(0,1),(1,1)]):
    assert poly==sum((bell[k]*mats[k][i][j]*G[k] for k in range(d+1)),P([]))
scale=Q(r['slow_polynomial_reference_scale_exact']);eta=Q(r['slow_polynomial_constant_subtraction_exact'])
assert scale>0 and eta>0 and polys[0]==scale*refpoly-eta
strict=Q(r['cap_margin_exact']);cut=Q(r['cap_cutoff_exact'])
def variation(poly,point):
    seq=[poly,poly.derivative()]
    while seq[-1].degree()>0:
        rem=-(seq[-2]%seq[-1]);assert rem
        z=rem.numer();seq.append(P(z/abs(z.content())))
    vals=[p(point) for p in seq];sg=[1 if z>0 else -1 for z in vals if z]
    return sum(x!=y for x,y in zip(sg,sg[1:]))
poly=polys[0]+polys[1]+strict
assert poly(-1)<0 and poly(cut)<0 and variation(poly,Q(-1))==variation(poly,cut)
Babs=[sum((bell[k]*abs(mats[k][i][j]) for k in range(d+1)),Q(0)) for i,j in [(0,0),(0,1),(1,1)]]
assert Babs==list(map(Q,r['angular_absolute_bounds']))
Tmax,Umax=Q(r['compact_T_upper_exact']),Q(r['compact_u_upper_exact'])
assert slow>=40 and gamma>0 and Umax>=2 and Tmax>1 and strict<eta
tailT=6*arb(Babs[1])*(-arb(gamma)*(Tmax-1)).exp()+3*arb(Babs[2])*(-2*arb(gamma)*(Tmax-1)).exp()
tailU=(2*arb(Babs[1]+Babs[2])*(-arb(gamma)*(Umax*Umax-1)).exp()
      +arb(Babs[0]+Babs[1])*(-arb(slow)*(2*cut*Umax-1)).exp()
      +arb(Babs[0]+2*Babs[1]+Babs[2])*(-arb(slow)*(Umax*Umax-1)).exp())
assert tailT<arb(eta)/2 and tailU<arb(strict)/2
AP=[arb_poly([arb(z) for z in p.coeffs()]) for p in polys]
cache={};evaluated=0
def translated(z):
    return [p(arb_poly([z.mid(),z.rad()])) for p in AP]
def geom(lo,hi,ul,uh):
    key=(lo,hi,ul,uh)
    if key in cache:return cache[key]
    T=arb(lo).union(arb(hi));U=arb(ul).union(arb(uh))
    wl=1+ul*ul-2*hi*ul;wh=1+uh*uh-2*lo*uh
    assert wh>=wl>=1
    W=arb(wl).union(arb(wh)).sqrt()
    products=[lo*ul,lo*uh,hi*ul,hi*uh]
    pl=max(Q(0),1-2*max(products));ph=1-2*min(products)
    out=(wl,wh,pl,ph,translated(T),translated(((1-U*T)/W).intersection(arb(0,1))),
         translated(((U-T)/W).intersection(arb(0,1))))
    if len(cache)<20000:cache[key]=out
    return out
def exp_interval(lo,hi):return arb(lo).exp().union(arb(hi).exp())
def combination(coeff,x,y):
    total=coeff[0]+(x+y)*coeff[1]+x*y*coeff[2]
    error=sum((abs(z) for z in total.coeffs()[1:]),arb(0))
    return total[0]+arb(0,error.upper())
def bound(box):
    global evaluated
    lo,hi,ul,uh,tl,th,dt,du,dT=box;evaluated+=1
    if 2*lo*ul>1:return None
    wl,wh,pl,ph,g0,g1,g2=geom(lo,hi,ul,uh)
    rA=exp_interval(-gamma*(th-1),-gamma*(tl-1))
    rB=exp_interval(-gamma*(th*uh*uh-1),-gamma*(tl*ul*ul-1))
    rC=exp_interval(-gamma*(th*wh-1),-gamma*(tl*wl-1))
    w1=exp_interval(-slow*th*ph,-slow*tl*pl)
    w2=exp_interval(-slow*th*(wh-1),-slow*tl*(wl-1))
    return combination(g0,rA,rB)+w1*combination(g1,rA,rC)+w2*combination(g2,rB,rC)
def actual(t,u,T):
    t,u,T=arb(t),arb(u),arb(T)
    w2=1+u*u-2*t*u;w=w2.sqrt()
    rA=(-arb(gamma)*(T-1)).exp();rB=(-arb(gamma)*(T*u*u-1)).exp();rC=(-arb(gamma)*(T*w2-1)).exp()
    def term(z,x,y):return AP[0](z)+(x+y)*AP[1](z)+x*y*AP[2](z)
    return term(t,rA,rB)+(-arb(slow)*T*(w2-u*u)).exp()*term((1-u*t)/w,rA,rC)+(-arb(slow)*T*(w2-1)).exp()*term((u-t)/w,rB,rC)
stack=[((Q(-1),Q(1,2),Q(1),Umax,Q(1),Tmax,0,0,0),None)]
passed=discarded=processed=maxdepth=0;start=time.time()
while stack:
    box,known=stack.pop();processed+=1
    lo,hi,ul,uh,tl,th,dt,du,dT=box
    if 2*lo*ul>1:discarded+=1;continue
    val=known if known is not None else bound(box)
    assert val is not None and val.is_finite()
    if val<0:
        passed+=1;maxdepth=max(maxdepth,dt,du,dT);continue
    if processed>args.budget or max(dt,du,dT)>=args.maxdepth:
        (D/'compact_inconclusive.json').write_text(json.dumps({'box':list(map(str,box)),'bound':str(val),'processed':processed},indent=2)+'\n')
        raise RuntimeError('Full compact proof subdivision budget exhausted')
    mt,mu,mT=(lo+hi)/2,(ul+uh)/2,(tl+th)/2
    if 2*mt*mu<=1:
        check=actual(mt,mu,mT)
        if check>0:
            point=[float(mt),float(mu),float(mT)]
            cutsfile=D/'cuts.json';cuts=json.loads(cutsfile.read_text()) if cutsfile.exists() else {}
            cuts['triangles']=cuts.get('triangles',[])+[point];cutsfile.write_text(json.dumps(cuts,indent=2)+'\n')
            (D/'compact_counterexample.json').write_text(json.dumps({'point_exact':list(map(str,[mt,mu,mT])),'value_interval':str(check)},indent=2)+'\n')
            raise RuntimeError('Positive compact triangle counterexample')
    choices=[
       [(lo,mt,ul,uh,tl,th,dt+1,du,dT),(mt,hi,ul,uh,tl,th,dt+1,du,dT)],
       [(lo,hi,ul,mu,tl,th,dt,du+1,dT),(lo,hi,mu,uh,tl,th,dt,du+1,dT)],
       [(lo,hi,ul,uh,tl,mT,dt,du,dT+1),(lo,hi,ul,uh,mT,th,dt,du,dT+1)]]
    scored=[]
    for axis,children in enumerate(choices):
        results=[bound(c) for c in children]
        top=max([float(v.upper()) for v in results if v is not None],default=-float('inf'))
        scored.append((top,[dt,du,dT][axis],children,results))
    _,_,children,results=min(scored,key=lambda z:(z[0],z[1]))
    stack.extend(zip(children,results))
    if processed%10000==0:print('COMPACT',processed,passed,discarded,'depth',maxdepth,'seconds',time.time()-start,flush=True)
out={'pass':True,'dimension':n,'nonorigin_blocks_PSD_exact':True,'origin_schur_identity_exact':True,
     'slow_reference_identity_exact':True,'mixed_cap_exact_Sturm':True,
     'large_T_tail_exact_interval':True,'large_u_tail_exact_interval':True,
     'complete_compact_domain':'-1<=t<=1/2,1<=u<=U,1<=T<=Tmax,2tu<=1',
     'compact_T_upper_exact':str(Tmax),'compact_u_upper_exact':str(Umax),
     'accepted_boxes':passed,'discarded_sorting_boxes':discarded,'processed':processed,
     'evaluated_bounds':evaluated,'max_axis_depth':maxdepth,'seconds':time.time()-start,
     'reference_Gaussian_angular_replay_required':True,
     'scope':'Matrix kernel only; radial Fourier/pair gates and density objective remain required'}
(D/'kernel_verified.json').write_text(json.dumps(out,indent=2)+'\n')
print('FULL MATRIX KERNEL PASS',json.dumps(out),flush=True)
