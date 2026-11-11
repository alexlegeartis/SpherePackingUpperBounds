"""Exact Delsarte inputs for midpoint and transverse lens projections."""
if not __debug__:
    raise RuntimeError('Do not disable assertions.')
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sys

from flint import fmpq as Q,fmpq_poly as Poly

HERE=Path(__file__).resolve().parent
R=F(27,25)
EDGES=list(map(F,['3/2','31/20','8/5','41/25']))


def projection(n,d,kind):
    upper=F(41,25)
    assert F(3,2)<=d<=upper
    U=R*R-d*d/4
    assert 0<U<1
    if kind=='midpoint':
        D=F(0);L=F(0);dim=n
        s=1-F(1,2)/U
        # Monotonicity holds for any two positive norms <=sqrt(U)<1.
        condition=1-U
        uniform_lower=F(0)
        remainder=F(0)
    else:
        assert kind=='transverse'
        D=(R*R-1)**2/(d*d)
        L=1-d*d/4-D/4
        assert L>0
        condition=L-U+1-D
        assert condition>0
        s=1-(1-D)/(2*U)
        dim=n-1
        A=(R*R-1)**2
        uniform_lower=1-upper*upper/4-A/(4*upper*upper)
        remainder=A*(4*R*R-A)
        assert uniform_lower>0 and remainder>0 and d**4>A
    return dict(packing_dimension=n,base_lower=str(d),projection=kind,
                ambient_dimension=dim,max_inner_product=str(s),
                transverse_longitudinal_difference_squared=str(D),
                norm_squared_lower=str(L),norm_squared_upper=str(U),
                monotonicity_margin=str(condition),valid_through_base=str(upper),
                uniform_norm_squared_lower=str(uniform_lower),
                angular_monotonicity_positive_remainder=str(remainder),
                zero_midpoint_case='If present, all other vectors have norm<1, so M<=1.')


def exact_certificate(dim,cs,s):
    X=Poly([0,1]);ps=[Poly([1]),X]
    for k in range(1,len(cs)-1):
        ps.append(((2*k+dim-2)*X*ps[-1]-k*ps[-2])/Q(k+dim-2))
    coeff=list(map(lambda x:Q(str(F(str(x)))),cs))
    assert coeff[0]==1 and all(v>=0 for v in coeff)
    P=sum((v*p for v,p in zip(coeff,ps)),Poly([]))
    sequence=[P,P.derivative()]
    while sequence[-1].degree()>0:
        remainder=-(sequence[-2]%sequence[-1])
        if not remainder:break
        z=remainder.numer()
        sequence.append(Poly(z/abs(z.content())))
    def variation(values):
        signs=[1 if x>0 else -1 for x in values if x]
        return sum(a!=b for a,b in zip(signs,signs[1:]))
    assert P(-1)<0 and P(s)<0
    assert variation([p(-1)for p in sequence])-variation([p(s)for p in sequence])==0
    bound=sum(coeff,Q(0));cap=int(bound)
    assert cap<=bound<cap+1
    return dict(gegenbauer_coefficients=list(map(str,coeff)),polynomial=list(map(str,P.coeffs())),
                rational_bound=str(bound),integer_bound=cap,exact_global_sign='passed')


def solve(dim,s,degree=20):
    import numpy as np
    from scipy.optimize import linprog
    from scipy.special import eval_gegenbauer,poch,comb,factorial
    lam=(dim-2)/2
    mats=[]
    for left,right in zip(np.linspace(-1,float(s),81)[:-1],np.linspace(-1,float(s),81)[1:]):
        jets=np.zeros((degree+1,degree+1));bc=np.zeros_like(jets)
        for k in range(degree+1):
            for ell in range(k,degree+1):
                jets[k,ell]=2**k*poch(lam,k)*eval_gegenbauer(ell-k,lam+k,left)/eval_gegenbauer(ell,lam,1)/factorial(k)
        for i in range(degree+1):
            for k in range(i+1):bc[i,k]=comb(i,k)/comb(degree,k)*(right-left)**k
        mats.append(bc@jets)
    M=np.vstack(mats)
    result=linprog(np.ones(degree),A_ub=M[:,1:],b_ub=-M[:,0]-1e-6,
                   bounds=[(0,None)]*degree,method='highs')
    assert result.success,result.message
    return exact_certificate(dim,[1.]+list(result.x),Q(str(s)))


def replay(path):
    data=json.loads(path.read_text())
    for row in data['certificates']:
        expected=projection(row['packing_dimension'],F(row['base_lower']),row['projection'])
        assert all(row[k]==v for k,v in expected.items())
        certificate=exact_certificate(row['ambient_dimension'],row['gegenbauer_coefficients'],Q(row['max_inner_product']))
        assert all(row[k]==v for k,v in certificate.items())
    return dict(verified=True,certificates=len(data['certificates']))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',type=Path)
    ap.add_argument('--output',type=Path,default=HERE/'spherical_lens_caps.json')
    args=ap.parse_args()
    if args.verify:
        out=replay(args.verify);print(out)
        args.verify.with_name(args.verify.stem+'_replay.json').write_text(json.dumps(out,indent=2)+'\n')
    else:
        rows=[]
        for n in [5,6,7]:
            for d in EDGES:
                for kind in ['midpoint','transverse']:
                    row=projection(n,d,kind)
                    row.update(solve(row['ambient_dimension'],F(row['max_inner_product'])))
                    rows.append(row)
                    print(n,d,kind,row['integer_bound'],float(Q(row['rational_bound'])),flush=True)
        args.output.write_text(json.dumps(dict(short_radius=str(R),certificates=rows),indent=2)+'\n')
