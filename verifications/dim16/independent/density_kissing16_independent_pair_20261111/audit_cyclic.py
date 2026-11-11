#!/usr/bin/env python3
"""Independent factor-form cyclic reconstruction and exact cache enclosure audit.

Reads the separately generated cyclic artifacts; writes only in this directory.
Does not itself certify global cyclic nonpositivity.
"""
if not __debug__:
    raise RuntimeError("Certificate verification requires assertions; do not run Python with -O or -OO.")

import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
import json
import hashlib
import math
from pathlib import Path
from fractions import Fraction
from time import monotonic
import numpy as np
from flint import fmpq,fmpq_mpoly_ctx,ctx
import verify
ctx.threads=1
HERE=Path(__file__).resolve().parent
CYCLIC=HERE.parent/'density_kissing16_independent_cyclic_20261111'


def main():
    start=monotonic()
    for name,digest in verify.EXPECTED.items():
        assert hashlib.sha256((verify.FROZEN/name).read_bytes()).hexdigest()==digest
    factors=json.loads((verify.FROZEN/'factors.json').read_text(),parse_float=str)
    candidate=json.loads((verify.FROZEN/'candidate.json').read_text(),parse_float=str)
    blocks,_=verify.psd_blocks(factors)
    th=list(map(verify.Q,candidate['record']['theta']))
    rr=list(map(verify.Q,candidate['record']['local_angular_coefficients']))
    mc=fmpq_mpoly_ctx.get(['u','v','t']);u,v,t=mc.gens()
    zero=mc.constant(0)
    def sub_univariate(p,x):
        return sum((c*x**i for i,c in enumerate(p.coeffs())),zero)
    def kernel(a,b,c):
        g=c-a*b;z=(1-a*a)*(1-b*b)
        qs=[mc.constant(1),g]
        for k in range(1,16):
            qs.append(((2*k+13)*g*qs[-1]-k*z*qs[-2])/fmpq(k+13))
        out=mc.constant(blocks[16][2][1][1])
        for k in range(17):
            if not th[k]:continue
            L,lam,_=blocks[17+k]
            p=verify.polys(16+2*k,16-k)
            pa=[sub_univariate(q,a) for q in p]
            pb=[sub_univariate(q,b) for q in p]
            # Direct bilinear factor expansion; no dense F or basis matrix.
            local=sum((sum((x*y for x,y in zip(row,pa)),zero)*sum((x*y for x,y in zip(row,pb)),zero) for row in L),zero)
            local+=lam*sum((x*y for x,y in zip(pa,pb)),zero)
            out+=th[k]*qs[k]*local
        return out
    poly=kernel(u,v,t)+kernel(u,t,v)+kernel(v,t,u)
    pglobal=verify.polys(16,24)
    R=sum((rr[j-1]*pglobal[j] for j in range(1,25)),verify.fmpq_poly([]))
    poly+=sub_univariate(R,u)+sub_univariate(R,v)+sub_univariate(R,t)
    raw=(CYCLIC/'polynomial.json').read_bytes()
    claimed=json.loads(raw)
    dct={tuple(e):verify.Q(c) for e,c in claimed['terms']}
    assert dct==poly.to_dict()
    assert poly==poly.compose(v,u,t)==poly.compose(u,t,v)
    assert all(claimed['source_sha256'][name]==digest for name,digest in verify.EXPECTED.items())
    print('Independent factor-form cyclic polynomial equality PASS:',len(dct),'terms',flush=True)
    degree=max(max(e) for e in dct);n=degree+1
    cap=verify.Q(candidate['cap_exact']);a=fmpq(-1);width=cap+1
    # Separate affine power substitution and monomial -> Bernstein conversion.
    affine=[[fmpq(math.comb(e,j))*a**(e-j)*width**j if j<=e else fmpq(0) for e in range(n)] for j in range(n)]
    bern=[[fmpq(math.comb(i,j),math.comb(degree,j)) if j<=i else fmpq(0) for j in range(n)] for i in range(n)]
    allindices=list(np.ndindex((n,n,n)))
    coeff={e:dct.get(e,fmpq(0)) for e in allindices}
    def transform(coeff, matrix, axis):
        out={}
        for e in allindices:
            idx=list(e)
            values=[]
            for j in range(n):
                if not matrix[e[axis]][j]:continue
                idx[axis]=j
                values.append(matrix[e[axis]][j]*coeff[tuple(idx)])
            out[e]=sum(values,fmpq(0))
        return out
    for axis in range(3):
        coeff=transform(coeff,affine,axis)
        coeff=transform(coeff,bern,axis)
    (HERE/'exact_cyclic_bernstein.json').write_text(json.dumps({'degree':degree,'coefficients':[str(coeff[e]) for e in allindices]},separators=(',',':'))+'\n')
    raw_npz=(CYCLIC/'bernstein_initial.npz').read_bytes()
    cache=np.load(CYCLIC/'bernstein_initial.npz')
    lo,hi=cache['lo'],cache['hi']
    assert lo.shape==hi.shape==(n,n,n)
    assert np.isfinite(lo).all() and np.isfinite(hi).all()
    for e in allindices:
        assert verify.Q(Fraction(float(lo[e])))<=coeff[e]<=verify.Q(Fraction(float(hi[e]))), e
    print('All exact rational Bernstein coefficients contained in saved binary64 intervals:',len(coeff),flush=True)
    # Audit the imported split function on the actual cache: compare every child
    # coefficient with exact rational de Casteljau for each of the three axes.
    import importlib.util
    spec=importlib.util.spec_from_file_location('cyclic_implementation',CYCLIC/'bernstein.py')
    implementation=importlib.util.module_from_spec(spec);spec.loader.exec_module(implementation)
    for axis in range(3):
        ll,lh,rl,rh=implementation.split(lo,hi,axis)
        other=[a for a in range(3) if a!=axis]
        for i in range(n):
            for j in range(n):
                ix=[0,0,0];ix[other[0]]=i;ix[other[1]]=j
                exact=[]
                for k in range(n):
                    ix[axis]=k;exact.append(coeff[tuple(ix)])
                lhs=[exact[0]];rhs=[exact[-1]]
                while len(exact)>1:
                    exact=[(x+y)/2 for x,y in zip(exact,exact[1:])]
                    lhs.append(exact[0]);rhs.append(exact[-1])
                rhs.reverse()
                for k in range(n):
                    ix[axis]=k;e=tuple(ix)
                    assert verify.Q(Fraction(float(ll[e])))<=lhs[k]<=verify.Q(Fraction(float(lh[e])))
                    assert verify.Q(Fraction(float(rl[e])))<=rhs[k]<=verify.Q(Fraction(float(rh[e])))
    out={
        'status':'PASS for audited construction and interval enclosures; global search completion is not checked here',
        'independent_factor_form_matches_exact_cyclic_polynomial':True,
        'exact_permutation_symmetry':True,
        'terms':len(dct),
        'degree_each_axis':degree,
        'exact_initial_bernstein_coefficients_checked':len(coeff),
        'all_initial_binary64_intervals_contain_exact_coefficients':True,
        'all_first_subdivision_intervals_contain_exact_child_coefficients_on_all_axes':True,
        'exact_cyclic_bernstein_sha256':hashlib.sha256((HERE/'exact_cyclic_bernstein.json').read_bytes()).hexdigest(),
        'polynomial_sha256':hashlib.sha256(raw).hexdigest(),
        'bernstein_initial_npz_sha256':hashlib.sha256(raw_npz).hexdigest(),
        'audited_script_sha256':hashlib.sha256((CYCLIC/'bernstein.py').read_bytes()).hexdigest(),
        'runtime_seconds':monotonic()-start,
    }
    (HERE/'cyclic_audit_result.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2),flush=True)

if __name__=='__main__':main()
