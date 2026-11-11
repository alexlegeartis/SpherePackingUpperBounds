"""Exact slow-Gaussian correction certificate for a polynomial-Gaussian pair."""
import math
from fractions import Fraction
from flint import fmpq as Q, fmpq_poly as Poly, arb, ctx


def sturm(p):
    seq=[p,p.derivative()]
    while seq[-1].degree()>0:
        q=-(seq[-2]%seq[-1])
        if not q: break
        z=q.numer();seq.append(Poly(z/abs(z.content())))
    return seq


def variation(values):
    s=[1 if v>0 else -1 for v in values if v]
    return sum(a!=b for a,b in zip(s,s[1:]))


def count(seq,lo,hi=None):
    vl=variation([p(lo) for p in seq])
    vh=variation([p[p.degree()] if hi is None else p(hi) for p in seq])
    return vl-vh


def upper_on_interval(p,lo,hi):
    m=arb((lo+hi)/2);h=arb((hi-lo)/2)
    co=[arb(0)]*(p.degree()+1)
    for z in reversed(p.coeffs()):
        for j in range(p.degree(),0,-1):co[j]=co[j]*m+co[j-1]
        co[0]=co[0]*m+arb(z)
    return co[0]+sum((abs(co[j])*h**j for j in range(1,len(co))),arb(0))


def verify_correction(r,P,H,r2,R2,pi_lo):
    assert r['kernel_family']=='solid_harmonic_cyclic_slow_gaussian'
    g=r['slow_gaussian_correction'];n=int(r['dimension']);a=Q(n,2)
    beta=Q(g['beta']);c=Q(g['c']);eta=Q(g['eta'])
    assert beta==Q(1,4) and eta>0 and arb(c)<arb.pi()*arb(r2)
    factor=Q(2**n) # beta^(-n/2), exact even for odd n
    negative_cutoff=a*beta-c*beta*beta
    assert negative_cutoff>0
    cost=eta*factor*(a/beta-c)
    shifted=H-cost;hs=sturm(shifted)
    assert shifted(0)>0 and shifted(negative_cutoff)>0
    assert count(hs,Q(0),negative_cutoff)==0
    # H itself is checked positive globally by the main verifier. Below the
    # cutoff the correction, after dividing by exp(-Y), is at least -cost.
    # Above it the correction is nonnegative.
    print('Repaired global Fourier positivity: PASS; origin loss',arb(cost),flush=True)
    intervals=[tuple(map(Q,x)) for x in g['repair_intervals']]
    assert intervals and all(lo<hi for lo,hi in intervals)
    start=R2*pi_lo;assert start>c
    ps=sturm(P);left=start
    for lo,hi in intervals:
        assert left<lo and lo>c
        assert P(left)<0 and P(lo)<0 and count(ps,left,lo)==0
        up=upper_on_interval(P,lo,hi)
        allowance=arb(eta)*arb(lo-c)*(arb(1-beta)*arb(lo)).exp()
        assert up<allowance,('REPAIR_INTERVAL_FAIL',str(lo),str(hi),str(up),str(allowance))
        print('Tail repair interval',lo,hi,'PASS; upper',up,'allowance',allowance,flush=True)
        left=hi
    assert P(left)<0 and count(ps,left,None)==0
    print('Repaired entire physical tail: PASS',flush=True)
    normalization=H[0]+eta*factor*(c-a/beta)
    assert normalization>0 and normalization==Q(g['corrected_normalization'])
    return eta,c,normalization
