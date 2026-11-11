"""Replay exact auxiliary tetrahedral inequalities; this is not a density verifier.

Only python-flint is needed. No optimization or floating-point sign decisions.
Usage: python verify_tetra_envelopes.py [certificate.json ...]
"""
import json
import sys
from pathlib import Path
from fractions import Fraction
from flint import fmpq as Q, fmpq_poly as Poly


def gegenbauer(d, degree):
    x=Poly([0,1]); out=[Poly([1])]
    if degree: out.append(x)
    for k in range(1, degree):
        out.append(((2*k+d-2)*x*out[-1]-k*out[-2])/(k+d-2))
    return out


def negative_on(poly, lo, hi):
    assert lo<hi and poly(lo)<0 and poly(hi)<0
    if poly.degree()==0: return
    chain=[poly,poly.derivative()]
    while chain[-1].degree()>0:
        rem=-(chain[-2]%chain[-1])
        if not rem: break
        chain.append(rem/abs(rem.leading_coefficient()))
    def variations(at):
        s=[1 if p(at)>0 else -1 for p in chain if p(at)]
        return sum(a!=b for a,b in zip(s,s[1:]))
    assert variations(lo)==variations(hi), 'Interior real zero blocks strict negativity'


def replay(path):
    data=json.loads(path.read_text()); assert data
    for row in data:
        assert row['geometry_version'] in ('sharp_common_neighbor_v1','conservative_projection_v1')
        sharp=row['geometry_version']=='sharp_common_neighbor_v1'
        assert row['envelope_family'] in ('affine','normal_harmonics')
        n=row['dimension'];q=Q(row['q'])
        assert isinstance(n,int) and n>=4 and 1<q<=Q(6,5)
        low_radius=1-q/4;upper_radius=q-Q(1,4)
        lo=1-q/(2*low_radius)
        delta=Q(0) if sharp else (q-1)**2
        hi=1-(1-delta)/(2*upper_radius);far=1-(q-delta)/(2*upper_radius)
        assert [Q(row[k])for k in ('lo','hi','far')]==[lo,hi,far]
        cs=list(map(Q,row['coefficients']));assert cs[0]==1 and all(c>=0 for c in cs)
        basis=gegenbauer(n-1,len(cs)-1)
        P=sum((c*p for c,p in zip(cs,basis)),Poly([]))
        assert P==Poly(list(map(Q,row['polynomial'])))
        B=sum(cs,Q(0));assert P(1)==B and B>0 and B==Q(row['B_exact'])
        if row['envelope_family']=='normal_harmonics':
            aa=list(map(Q,row['envelope_coefficients']));assert all(a>=0 for a in aa)
            leg=gegenbauer(3,len(aa))
            E=sum((a*((-1)**(ell+1)*leg[ell]-Q(1,3))for ell,a in enumerate(aa,1)),Poly([]))
            assert E==Poly(list(map(Q,row['envelope_polynomial'])))
        else:
            a=Q(row['slope']);assert a>=0
            E=a*Poly([Q(-1,3),1])
        negative_on(P,Q(-1),far);negative_on(P-E,lo,hi)
        f=Fraction(str(B));den=10**8;ceil=(f.numerator*den+f.denominator-1)//f.denominator
        upper=f'{ceil//den}.{ceil%den:08d}'
        print(f'PASS n={n}, q={q}: aggregate 3 T <= B E, B <={upper}; {path.name}',flush=True)
    return len(data)

if __name__=='__main__':
    root=Path(__file__).resolve().parent
    paths=list(map(Path,sys.argv[1:])) or [root/'affine_baselines.json']
    count=sum(replay(p)for p in paths)
    print(f'ALL {count} exact auxiliary certificates PASS. No density claim.')
