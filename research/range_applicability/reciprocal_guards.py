"""Exact reciprocal-Gaussian guards for physical and Fourier tail islands.

The physical additions are eta(c-X)e^(-X/4) and zeta(c-X)e^(-4X).
This module never treats small sign violations as harmless.
"""
from flint import fmpq as Q,fmpq_poly as Poly,arb
from slow_gaussian import sturm,count
from math import comb


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
    assert r['kernel_family']=='solid_harmonic_cyclic_reciprocal_guards'
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
