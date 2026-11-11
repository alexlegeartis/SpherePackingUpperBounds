"""Kernel-independent exact reciprocal-Gaussian guards.

Adapted from the repository's solid-cyclic verifier. The analytic guard
identities apply to any polynomial-Gaussian radial pair; no kernel sign
condition is inferred from these checks.
"""
if not __debug__: raise RuntimeError("Verification must not run under python -O.")

import math
from math import comb
from fractions import Fraction
from flint import fmpq as Q, fmpq_poly as Poly, arb
PI_LO=Q('314159265358979323846264338327950288419716939937510/100000000000000000000000000000000000000000000000000')

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



def bernstein_upper(p,lo,hi):
    d=p.degree();shift=p(Poly([lo,1]));width=hi-lo
    return max(sum((shift[j]*width**j*Q(comb(i,j),comb(d,j)) for j in range(i+1)),Q(0)) for i in range(d+1))


def certify_complement(p,start,intervals,sign):
    seq=sturm(p);left=start
    for lo,hi in intervals:
        assert start<=lo<hi and left<=lo
        if left<lo:
            assert sign*p(left)>0 and sign*p(lo)>0 and count(seq,left,lo)==0
        left=hi
    assert sign*p(left)>0 and count(seq,left,None)==0


def verify_guards(r,P,H,r2,R2,pi_lo):
    assert r['gaussian_guard_family']=='reciprocal_gaussians'
    n=int(r['dimension']);g=r['reciprocal_gaussian_guards']
    c=Q(g['c']);eta=Q(g['eta']);zeta=Q(g['zeta'])
    assert eta>=0 and zeta>=0 and c>Q(n,8) and arb(c)<arb.pi()*arb(r2)
    D=Q(2**n)*(2*n-c);E=(c-Q(n,8))/Q(2**n);y0=(2*n-c)/16
    assert D>0 and E>0 and y0>0
    # Exact differentiated Gaussian identities are stated in the proof note.
    # The slow Fourier guard can be negative only on [0,y0], by at most eta D.
    low=H-eta*D;seq=sturm(low)
    assert low(0)>0 and low(y0)>0 and count(seq,Q(0),y0)==0
    physical=[tuple(map(Q,z)) for z in g['physical_intervals']]
    fourier=[tuple(map(Q,z)) for z in g['fourier_intervals']]
    start=R2*pi_lo;assert start>c
    certify_complement(P,start,physical,-1)
    for lo,hi in physical:
        up=bernstein_upper(P,lo,hi)
        allowance=arb(eta)*arb(lo-c)*(arb(3)*arb(lo)/4).exp()
        assert arb(up)<allowance,('PHYSICAL_GUARD',str(lo),str(hi),str(arb(up)),str(allowance))
    # Beyond y0 the slow guard is nonnegative; the fast Fourier guard is
    # everywhere positive because c>n/8. Thus it suffices to fix -H islands.
    certify_complement(H,y0,fourier,1)
    for lo,hi in fourier:
        assert lo>=y0
        up=bernstein_upper(-H,lo,hi)
        allowance=arb(zeta/Q(2**n))*arb(c-Q(n,8)+lo/16)*(arb(3)*arb(lo)/4).exp()
        assert arb(up)<allowance,('FOURIER_GUARD',str(lo),str(hi),str(arb(up)),str(allowance))
    normalization=H[0]-eta*D+zeta*E
    assert normalization>0 and normalization==Q(g['corrected_normalization'])
    print('Reciprocal Gaussian identities and both global signs: PASS',len(physical),len(fourier),flush=True)
    return c*(eta+zeta),normalization

def islands(p,start,desired_sign):
    real=[]
    for z,m in p.complex_roots():
        if z.imag==0 and float(z.real)>float(start):real.append(float(z.real))
    real.sort();points=[float(start)]+real;bad=[]
    for i in range(len(points)-1):
        mid=Q(str(Fraction(str((points[i]+points[i+1])/2))))
        if desired_sign*p(mid)<0:
            lo=start if i==0 else Q(math.floor(points[i]*1000)-1,1000)
            hi=Q(math.ceil(points[i+1]*1000)+1,1000)
            bad.append((max(start,lo),hi))
    merged=[]
    for lo,hi in bad:
        if merged and lo<=merged[-1][1]:merged[-1]=(merged[-1][0],max(hi,merged[-1][1]))
        else:merged.append((lo,hi))
    certify_complement(p,start,merged,desired_sign)
    return merged


def outward(x):
    if x<=0:return Q(0)
    value=float(x.upper());exp=math.floor(math.log10(value))-1
    step=Q(10)**exp;q=Q(math.ceil(value/float(step)))*step
    while not arb(q)>x:q+=step
    return q


