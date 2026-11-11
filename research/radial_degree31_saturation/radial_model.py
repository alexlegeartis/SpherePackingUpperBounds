"""Exact model data and validated genuine sample rows for radial duals."""
import json
import math
from fractions import Fraction
from pathlib import Path
from flint import arb, arb_poly, fmpq as Q, fmpq_poly as P

if not __debug__:
    raise RuntimeError('Proof checks require assertions; do not use Python -O')
D=Path(__file__).resolve().parent


class Model:
    def __init__(self,n,B):
        self.n=n
        self.B=Q(B)
        assert n in [46,47] and self.B>0
        raw=json.loads((D/'inputs/source48_radial.json').read_text())
        normalizer=Q(raw['fourier_polynomial'][0])
        assert normalizer>0
        self.source=P(list(map(Q,raw['physical_polynomial'])))/normalizer
        alpha=Q(n,2)
        x=P([0,1])
        L0,L1=P([1]),alpha-x
        monomial_transforms=[L0,L1]
        self.fourier=self.source[0]*L0+self.source[1]*L1
        fact=1
        for k in range(1,self.source.degree()):
            L2=((2*k+alpha-x)*L1-(k+alpha-1)*L0)/(k+1)
            fact*=k+1
            monomial_transforms.append(fact*L2)
            self.fourier+=self.source[k+1]*fact*L2
            L0,L1=L1,L2
        self.H0=self.fourier(0)
        assert self.H0>0
        self.P0=self.source(0)
        self.c=[self.P0-self.H0]+[Q(1-(-1)**k) for k in range(1,32)]+[Q(1)]
        self.AP=arb_poly([arb(c) for c in self.source.coeffs()])
        self.AH=arb_poly([arb(c) for c in self.fourier.coeffs()])
        self.alpha=alpha
        # Verify all 32 normalized Laguerre Fourier eigenvalue identities
        # in exact rational polynomial arithmetic, including odd dimensions.
        E0,E1=P([1]),alpha-2*x
        eigen=[E0,E1/E1(0)]
        for k in range(1,31):
            E2=((2*k+alpha-2*x)*E1-(k+alpha-1)*E0)/(k+1)
            eigen.append(E2/E2(0))
            E0,E1=E1,E2
        for k,e in enumerate(eigen):
            transform=sum((e[j]*monomial_transforms[j] for j in range(k+1)),P([]))
            assert transform==(-1)**k*e and e(0)==1

    def lag(self,X):
        L0,L1=arb(1),arb(self.alpha)-2*X
        out=[L0,L1/arb(self.alpha)]
        normalizer=self.alpha
        for k in range(1,31):
            L2=((2*k+self.alpha-2*X)*L1-(k+self.alpha-1)*L0)/(k+1)
            normalizer*=Q(k+self.alpha,k+1)
            out.append(L2/arb(normalizer))
            L0,L1=L1,L2
        return out

    def row(self,spec,scaled=True):
        kind=spec['kind']
        assert kind in ['pair','physical','fourier']
        if kind=='pair':
            u=Q(spec['u'])
            assert 6<=u<=7
            X=arb.pi()*arb(u)
        elif spec.get('pi_interpolation') is not None:
            fraction=Q(spec['pi_interpolation'])
            assert 0<=fraction<=1 and kind=='physical'
            X=7*arb.pi()+(200-7*arb.pi())*arb(fraction)
        else:
            X=arb(Q(spec['X']))
            if kind=='physical':assert X>=7*arb.pi()
            if kind=='fourier':assert X>=0
        lag=self.lag(X)
        if kind=='fourier':
            row=[-self.AH(X)]+[-(-1)**k*lag[k] for k in range(32)]+[arb(0)]
        else:
            last=arb(0)
            if kind=='pair':
                h=(-arb(Q(17,8))*arb.pi()*arb(u-6)).exp()
                last=-X.exp()*(2*h-h*h)/arb(self.B)
            row=[self.AP(X)]+lag+[last]
        if scaled:
            multiplier=Q(spec['multiplier_exact'])
            assert multiplier>0
            row=[arb(multiplier)*x for x in row]
        A0=row[1]
        transformed=[row[0]-A0*arb(self.H0)]+[row[k+1]-A0*(-1)**k for k in range(1,32)]+[row[-1]]
        return transformed,A0,row


def rational_float(x):
    return str(Fraction(float(x)).limit_denominator(10**14))


def sample_specs(n):
    import numpy as np
    specs=[]
    for i in range(1001):specs.append({'kind':'pair','u':str(Q(6)+Q(i,1000))})
    for i in range(2001):specs.append({'kind':'physical','pi_interpolation':str(Q(i,2000))})
    for X in np.geomspace(200,2000,401):specs.append({'kind':'physical','X':rational_float(X)})
    for i in range(2401):specs.append({'kind':'fourier','X':str(Q(i,12))})
    for X in np.geomspace(200,2000,401):specs.append({'kind':'fourier','X':rational_float(X)})
    extra=json.loads((D/f'inputs/extra{n}.json').read_text())
    for kind,points in extra.items():
        assert kind in ['physical','fourier']
        for X in points:
            exact=Q(rational_float(X))
            if kind=='physical' and not arb(exact)>=7*arb.pi():continue
            specs.append({'kind':kind,'X':str(exact)})
    return specs
