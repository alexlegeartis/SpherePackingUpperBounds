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
AP=arb_poly([arb(z) for z in poly.coeffs()]);unit=arb(0).union(arb(1))
def angular_box(z):
    z=z.intersection(arb(0,1))
    terms=AP(arb_poly([z.mid(),z.rad()]))
    error=sum((abs(v) for v in terms.coeffs()[1:]),arb(0))
    return terms[0]+arb(0,error.upper())
def point_ratio(g,v):
    ratio=1+arb(g)/arb(v)
    power=ratio**int(-a) if a.denominator==1 else ratio**arb(-a)
    return power*(-arb(b)*arb(g)).exp()
def point_h(v):
    v=arb(v)
    power=v**int(-a) if a.denominator==1 else v**arb(-a)
    return power*(-arb(b)*(v-1)).exp()
stack=[(Q(-1),Q(1,2),Q(1),U,0,0)];passed=discarded=processed=maxdepth=0;start=time.time()
while stack:
    lo,hi,ul,uh,dt,du=stack.pop();processed+=1
    if 2*lo*ul>1:discarded+=1;continue
    t=arb(lo).union(arb(hi));u=arb(ul).union(arb(uh))
    wl,wh=1+ul*ul-2*hi*ul,1+uh*uh-2*lo*uh
    assert wh>=wl>=1
    w=arb(wl).union(arb(wh)).sqrt()
    products=[lo*ul,lo*uh,hi*ul,hi*uh]
    gl=max(Q(0),1-2*max(products));gh=1-2*min(products)
    # R(g,v)=(1+g/v)^(-a)exp(-bg) decreases in g. Its v derivative has sign a.
    vlow,vhigh=(uh*uh,ul*ul) if a<0 else (ul*ul,uh*uh)
    weight1=point_ratio(gh,vlow).union(point_ratio(gl,vhigh)).intersection(unit)
    weight2=point_h(wh).union(point_h(wl)).intersection(unit)
    first=angular_box(t)
    first_upper=min(first.upper(),(-arb(delta)).upper())
    bound=arb(first_upper)+weight1*angular_box((1-u*t)/w)+weight2*angular_box((u-t)/w)
    assert bound.is_finite()
    if bound<0:passed+=1;maxdepth=max(maxdepth,dt,du);continue
    assert max(dt,du)<35,(str(lo),str(hi),str(ul),str(uh))
    if dt<=du:
        mid=(lo+hi)/2;stack.extend([(lo,mid,ul,uh,dt+1,du),(mid,hi,ul,uh,dt+1,du)])
    else:
        mid=(ul+uh)/2;stack.extend([(lo,hi,ul,mid,dt,du+1),(lo,hi,mid,uh,dt,du+1)])
    if processed%100000==0:print('INDEPENDENT CYCLIC',processed,passed,'seconds',time.time()-start,flush=True)
out={'pass':True,'dimension':n,'degree':d,'precision_bits':ctx.prec,
     'closed_form_harmonics_exact':True,'constant_spherical_moment_exact':True,
     'exact_Bernstein_cap_leaves':cap_leaves,'exact_Bernstein_cap_depth':cap_depth,
     'unique_sign_transition_exact_Sturm':True,'complete_compact_domain_t':'[-1,1/2]',
     'compact_u_upper_exact':str(U),'certified_rectangles':passed,'discarded_sorting_rectangles':discarded,
     'processed':processed,'max_axis_depth':maxdepth,'monotone_ratio_endpoint_enclosures':True,
     'infinite_tail_rigorous':True,'all_scale_hypotheses_checked':True,
     'profile_a_exact':str(a),'profile_b_exact':str(b),'seconds':time.time()-start,
     'scope':'Independent angular/all-scale cyclic proof; independent radial proof required'}
(D/'independent_angular_result.json').write_text(json.dumps(out,indent=2)+'\n')
print('INDEPENDENT ANGULAR PASS',n,passed,maxdepth,flush=True)
