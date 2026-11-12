"""Exact-cone parametrization of g(v)=v phi'(v), v>=1.

Coordinates: g(1), jumps at knots, slopes on intervals [1,k0],...,[klast,infinity).
All coordinates nonnegative. Mathematical discovery support, not density proof.
"""
import numpy as np

def ramp_integral(v,a):
    """Integral_1^v (t-a)_+/t dt."""
    z=np.maximum(np.asarray(v)/a,1)
    return a*(z-1-np.log(z))

def basis(v, knots):
    v=np.asarray(v); k=np.asarray(knots)
    assert np.all(k>1) and np.all(np.diff(k)>0) and np.all(v>=1)
    ramps=[ramp_integral(v,a) for a in np.r_[1,k]]
    slopes=[ramps[i]-ramps[i+1] for i in range(len(k))]+[ramps[-1]]
    return np.stack([np.log(v)]+[np.maximum(np.log(v/a),0) for a in k]+slopes,axis=-1)

def old_to_new(p,VK,LK):
    """Convert signed linear/log hinges into monotone coordinates."""
    p=np.asarray(p);VK=np.asarray(VK);LK=np.asarray(LK);k=np.unique(np.r_[VK,LK]);nv=len(VK)
    a=p[0];b=p[1:1+nv];c=p[1+nv];d=p[2+nv:]
    jumps=np.array([t*b[VK==t].sum()+d[LK==t].sum() for t in k])
    slopes=np.array([a+b[VK<=t].sum() for t in np.r_[1,k]])
    return k,np.r_[a+c,jumps,slopes]

def new_to_old(q,k):
    """Aligned knots VK=LK=k. Returns signed hinge coefficients."""
    k=np.asarray(k);q=np.asarray(q);m=len(k);g0=q[0];J=q[1:1+m];A=q[1+m:]
    assert len(A)==m+1
    a=A[0];b=np.diff(A);c=g0-a;d=J-k*b
    return np.r_[a,b,c,d]

if __name__=='__main__':
    from scipy.integrate import quad
    rng=np.random.default_rng(9);k=np.array([1.003,1.006,1.05,1.2,1.4,2.])
    q=rng.uniform(0,10,2*len(k)+2);p=new_to_old(q,k);kk,qq=old_to_new(p,k,k)
    assert np.allclose(q,qq)
    v=np.r_[np.linspace(1,3,2001),np.geomspace(3,100,100)]
    nv=len(k);old=np.concatenate([(v-1)[:,None],np.maximum(v[:,None]-k,0),np.log(v)[:,None],np.maximum(np.log(v[:,None]/k),0)],axis=1)@p
    assert np.max(abs(old-basis(v,k)@q))<1e-10
    def g(t):
        return q[0]+sum(q[1+i] for i,x in enumerate(k) if t>x)+sum(q[1+len(k)+i]*max(0,min(t,b)-a) for i,(a,b) in enumerate(zip(np.r_[1,k],np.r_[k,np.inf])))
    for x in [1.02,1.7,5]:
        val=quad(lambda t:g(t)/t,1,x,points=list(k[k<x]),epsabs=1e-10)[0]
        assert abs(val-basis(np.array(x),k)@q)<1e-9
    print('PASS: conversion identity and quadrature checks (numerical tests, not density certificate)')

def knot_jacobian(v,q,k):
    """d phi/d k for fixed monotone coordinates, away from knot equality."""
    v=np.asarray(v); k=np.asarray(k);q=np.asarray(q);m=len(k);J=q[1:1+m];A=q[1+m:]
    return (A[:-1]-A[1:])*np.maximum(np.log(v[...,None]/k),0)-J/k*(v[...,None]>k)
