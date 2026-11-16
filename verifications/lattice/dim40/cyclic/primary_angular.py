"""Exact angular PSD/Sturm checks and full-domain log-Gaussian cyclic proof."""
import argparse,json,math,sys,time
from pathlib import Path
from flint import fmpq as Q,fmpq_poly as P,arb,arb_poly,ctx
from core import rational,normalized_gegenbauer
if not __debug__:raise RuntimeError('Assertions required')
sys.set_int_max_str_digits(0);ctx.prec=384
ap=argparse.ArgumentParser();ap.add_argument('package');ap.add_argument('--depth',type=int,default=34)
args=ap.parse_args();D=Path(args.package).resolve()
r=json.loads((D/'angular.json').read_text());cert=json.loads((D/'certificate.json').read_text())
n,d=r['dimension'],r['degree']
assert n==cert['dimension'] and rational(r['profile_a_exact'])==rational(cert['a_exact']) and rational(r['profile_b_exact'])==rational(cert['b_exact'])
a,b=rational(r['profile_a_exact']),rational(r['profile_b_exact'])
assert n>=3 and a+b>0 and b>0
poly=P(list(map(Q,r['polynomial_coefficients'])));cs=list(map(Q,r['gegenbauer_coefficients']))
G=normalized_gegenbauer(n,d);B=Q(r['B_exact']);cut=Q(r['cutoff']);delta=Q(r['negative_cap_delta'])
assert Q(1,2)<cut<1 and delta>0 and d==poly.degree() and B>0
assert len(cs)==d+1 and cs[0]==1 and min(cs)>=0 and sum(cs)==B==poly(1)
assert sum((c*g for c,g in zip(cs,G)),P([]))==poly
def sturm(p):
    seq=[p,p.derivative()]
    while seq[-1].degree()>0:
        rem=-(seq[-2]%seq[-1]);assert rem
        z=rem.numer();seq.append(P(z/abs(z.content())))
    return seq
def variation(seq,x):
    vals=[p(x) for p in seq];sg=[1 if z>0 else -1 for z in vals if z]
    return sum(x!=y for x,y in zip(sg,sg[1:]))
seq=sturm(poly)
assert poly(-1)<0 and poly(cut)<0 and poly(1)>0
assert variation(seq,Q(-1))==variation(seq,cut)
assert variation(seq,cut)-variation(seq,Q(1))==1
cap=poly+delta;seq=sturm(cap)
assert cap(-1)<0 and cap(Q(1,2))<0
assert variation(seq,Q(-1))==variation(seq,Q(1,2))
tail_u=int(r['tail_u'])
power_log=max(Q(0),-a)
def tail_bound(U):
    assert U*U>power_log/b
    return arb(B)*(arb(4)**arb(power_log)*(-arb(b)*(2*cut*U-1)).exp()+
                   arb(U)**arb(2*power_log)*(-arb(b)*(U*U-1)).exp())
while not tail_bound(tail_u)<arb(delta)/2:
    tail_u*=2
assert tail_u>=2 and b*(tail_u**2-1)>=b*(2*cut*tail_u-1)>0
# A single positive Taylor term gives a rational lower bound on exp.
z=b*(2*cut*tail_u-1);power=max(1,int(float(z)))
if power_log.denominator==1:
    p=int(power_log)
    assert z**power/math.factorial(power)>4*B*(4**p+Q(tail_u)**(2*p))/delta
print('EXACT PSD/CAP/TRANSITION/TAIL PASS',n,str(a),str(b),'tail_u',tail_u,flush=True)
AP=arb_poly([arb(c) for c in poly.coeffs()]);cache={}
def first_upper(lo,hi):
    if (lo,hi) not in cache:
        q=poly(P([lo,hi-lo]))
        bs=[sum((q[j]*Q(math.comb(i,j),math.comb(d,j)) for j in range(i+1)),Q(0)) for i in range(d+1)]
        cache[lo,hi]=max(bs)
    return arb(cache[lo,hi])
def angular_interval(z):
    z=z.intersection(arb(0,1))
    terms=AP(arb_poly([z.mid(),z.rad()]))
    err=sum((abs(c) for c in terms.coeffs()[1:]),arb(0))
    return terms[0]+arb(0,err.upper())
def nonnegative(z):
    return arb(0).union(z.upper()) if z.lower()<0 else z
stack=[(Q(-1),Q(1,2),Q(1),Q(tail_u),0,0)]
passed=discarded=processed=maxdepth=0;start=time.time()
while stack:
    lo,hi,ul,uh,dt,du=stack.pop();processed+=1
    if 2*lo*ul>1:discarded+=1;continue
    T=arb(lo).union(arb(hi));U=arb(ul).union(arb(uh))
    wlo=1+ul*ul-2*hi*ul;whi=1+uh*uh-2*lo*uh
    assert whi>=wlo>=1
    W2=arb(wlo).union(arb(whi));W=W2.sqrt()
    U2=arb(ul*ul).union(arb(uh*uh))
    gap1=nonnegative(1-2*T*U);gap2=nonnegative(U2-2*T*U)
    ratio1=1+1/U2-2*T/U
    if ratio1.lower()<1:ratio1=arb(1).union(ratio1.upper())
    weight1=(-arb(a)*ratio1.log()-arb(b)*gap1).exp()
    weight2=(-arb(a)*W2.log()-arb(b)*gap2).exp()
    # h decreases on v>=1 because a+b>0,b>0; both are allowed-radius ratios.
    unit=arb(0).union(arb(1))
    weight1=weight1.intersection(unit);weight2=weight2.intersection(unit)
    bound=first_upper(lo,hi)+weight1*angular_interval((1-U*T)/W)+weight2*angular_interval((U-T)/W)
    if bound.is_finite() and bound<0:
        passed+=1;maxdepth=max(maxdepth,dt,du);continue
    if max(dt,du)>=args.depth:
        point={'t_exact':str((lo+hi)/2),'u_exact':str((ul+uh)/2),
               't_box':[str(lo),str(hi)],'u_box':[str(ul),str(uh)],'bound':str(bound),'inconclusive':True}
        (D/'cyclic_failure.json').write_text(json.dumps(point,indent=2)+'\n')
        raise RuntimeError('Cyclic subdivision inconclusive')
    if first_upper(lo,hi)>=0 or dt<=du:
        mid=(lo+hi)/2;stack.extend([(lo,mid,ul,uh,dt+1,du),(mid,hi,ul,uh,dt+1,du)])
    else:
        mid=(ul+uh)/2;stack.extend([(lo,hi,ul,mid,dt,du+1),(lo,hi,mid,uh,dt,du+1)])
    if processed%25000==0:print('CYCLIC',processed,passed,discarded,'seconds',time.time()-start,flush=True)
out={'pass':True,'dimension':n,'degree':d,'profile_a_exact':str(a),'profile_b_exact':str(b),
     'tail_u':tail_u,'strict_cap_delta':str(delta),'complete_domain_t':'[-1,1/2]',
     'certified_rectangles':passed,'discarded_sorting_rectangles':discarded,'max_depth':maxdepth,
     'processed':processed,'tail_exponential_rational_lower_check':True,
     'all_scale_hypothesis':'a+b>0,b>0 imply v*phi_prime(v)=a+b*v is nonnegative and nondecreasing on v>=1',
     'seconds':time.time()-start,'scope':'Angular and cyclic only; radial signs/coupling and density normalization required'}
(D/'primary_angular_result.json').write_text(json.dumps(out,indent=2)+'\n')
print('FULL ANGULAR PASS',json.dumps(out),flush=True)
