#!/usr/bin/env python3
"""Independent fixed-point INTEGER proof replay of the completed Gram covering.

Floating-point arithmetic is used only to reproduce the discovered splitting
choices and lambda proposals. Every accepted or Gram-excluded leaf is checked
using integer upper enclosures on a fixed 2^-80 grid. Symmetry exclusion and
coverage use exact rational boxes. No floating-point sign can accept a leaf.
"""
if not __debug__:
    raise RuntimeError("Certificate verification requires assertions; do not run Python with -O or -OO.")

import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from time import monotonic
import numpy as np
from flint import arb, fmpq, ctx
ctx.prec=192;ctx.threads=1
HERE=Path(__file__).resolve().parent
CYCLIC=HERE.parent/'density_kissing16_independent_cyclic_20261111'
BITS=80
SCALE=1<<BITS


def ceil_scaled(q):
    q=Fraction(q)
    return -((-q.numerator*SCALE)//q.denominator)

def split_float_interval(lower,upper,axis):
    work_lower=np.moveaxis(lower,axis,0).copy()
    work_upper=np.moveaxis(upper,axis,0).copy()
    n=work_lower.shape[0]
    left_lower=[work_lower[0].copy()];left_upper=[work_upper[0].copy()]
    right_lower=[work_lower[-1].copy()];right_upper=[work_upper[-1].copy()]
    for unused in range(1,n):
        work_lower=np.nextafter(np.nextafter(work_lower[:-1]+work_lower[1:],-np.inf)/2,-np.inf)
        work_upper=np.nextafter(np.nextafter(work_upper[:-1]+work_upper[1:],np.inf)/2,np.inf)
        left_lower.append(work_lower[0].copy());left_upper.append(work_upper[0].copy())
        right_lower.append(work_lower[-1].copy());right_upper.append(work_upper[-1].copy())
    return tuple(np.moveaxis(np.stack(x),0,axis) for x in (left_lower,left_upper,list(reversed(right_lower)),list(reversed(right_upper))))

def split_integer_upper(upper,axis):
    work=np.moveaxis(upper,axis,0).copy()
    n=work.shape[0]
    left=[work[0].copy()];right=[work[-1].copy()]
    for unused in range(1,n):
        work=(work[:-1]+work[1:]+1)//2  # Exact ceil((a+b)/2), including a+b<0.
        left.append(work[0].copy());right.append(work[-1].copy())
    return np.moveaxis(np.stack(left),0,axis),np.moveaxis(np.stack(list(reversed(right))),0,axis)

def propose_multiplier(h,g):
    # Only reproduces a candidate. Exact integer acceptance is separate below.
    neg=g<0;pos=g>0;zero=~(neg|pos)
    if np.any(h[zero]>=0):return None
    lower=max(0.,float(np.max(h[neg]/(-g[neg])))) if np.any(neg) else 0.
    upper=float(np.min(-h[pos]/g[pos])) if np.any(pos) else np.inf
    if not(lower<upper and np.isfinite(lower)):return None
    lam=(lower+upper)*.5 if np.isfinite(upper) else np.nextafter(lower,np.inf)
    if not(lam>=0 and np.isfinite(lam)):return None
    combined=np.nextafter(h+np.nextafter(lam*g,np.inf),np.inf)
    return lam if np.max(combined)<0 else None


def main():
    start=monotonic()
    audit=json.loads((HERE/'cyclic_audit_result.json').read_text())
    for name,key in [('polynomial.json','polynomial_sha256'),('bernstein_initial.npz','bernstein_initial_npz_sha256'),('bernstein.py','audited_script_sha256')]:
        assert hashlib.sha256((CYCLIC/name).read_bytes()).hexdigest()==audit[key]
    expected=json.loads((CYCLIC/'bernstein_gram_result.json').read_text())
    assert expected['verified']
    exact_bytes=(HERE/'exact_cyclic_bernstein.json').read_bytes()
    assert hashlib.sha256(exact_bytes).hexdigest()==audit['exact_cyclic_bernstein_sha256']
    exact=json.loads(exact_bytes)
    degree=exact['degree'];n=degree+1;shape=(n,n,n)
    HI=np.array([ceil_scaled(q) for q in exact['coefficients']],dtype=object).reshape(shape)
    cached=np.load(CYCLIC/'bernstein_initial.npz');lo,hi=cached['lo'],cached['hi']
    for e,c in zip(np.ndindex(shape),exact['coefficients']):
        assert Fraction(c)<=Fraction(int(HI[e]),SCALE)
    # Independently form G as monomials and convert each coordinate to degree-d
    # Bernstein coefficients, retaining exact rationals and integer upper bounds.
    width=Fraction(15003,10000)
    P=[-1+width*Fraction(i,degree) for i in range(n)]
    Q=[1-2*width*Fraction(i,degree)+width**2*Fraction(i*(i-1),degree*(degree-1)) for i in range(n)]
    GI=np.empty(shape,dtype=object);gl=np.empty(shape);gh=np.empty(shape)
    for e in np.ndindex(shape):
        i,j,k=e
        value=1+2*P[i]*P[j]*P[k]-Q[i]-Q[j]-Q[k]
        GI[e]=ceil_scaled(value)
        ball=arb(fmpq(value.numerator,value.denominator))
        gl[e]=np.nextafter(float(ball.lower()),-np.inf)
        gh[e]=np.nextafter(float(ball.upper()),np.inf)
        assert Fraction(float(gl[e]))<=value<=Fraction(float(gh[e]))
        assert value<=Fraction(int(GI[e]),SCALE)
    stack=[(lo,hi,gl,gh,HI,GI,((Fraction(0),Fraction(1)),)*3,0)]
    counts={'checked':0,'passed_direct':0,'passed_gram_multiplier':0,'excluded_symmetry':0,'excluded_gram':0}
    maxdepth=0;covered=Fraction(0);transcript=hashlib.sha256();leaves=[];minmargin=None
    while stack:
        l,h,gl,gh,HI,GI,box,depth=stack.pop();counts['checked']+=1;maxdepth=max(maxdepth,depth)
        kind=None;lam=None
        if np.max(h)<0:
            # This EXACT integer sign, not the floating-point proposal, proves H<0.
            upper=int(np.max(HI));assert upper<0
            margin=Fraction(-upper,SCALE);minmargin=margin if minmargin is None else min(minmargin,margin)
            counts['passed_direct']+=1;kind='P'
        elif box[0][0]>box[1][1] or box[1][0]>box[2][1]:
            counts['excluded_symmetry']+=1;kind='S'
        elif np.max(gh)<0:
            assert int(np.max(GI))<0  # Exact proof G<0 throughout the box.
            counts['excluded_gram']+=1;kind='G'
        else:
            lam=propose_multiplier(h,gh)
            if lam is not None:
                num,den=lam.as_integer_ratio();assert num>=0 and den>0
                upper=int(np.max(HI*den+GI*num));assert upper<0
                margin=Fraction(-upper,SCALE*den);minmargin=margin if minmargin is None else min(minmargin,margin)
                counts['passed_gram_multiplier']+=1;kind='M'
        if kind is not None:
            covered+=math.prod(b-a for a,b in box)
            transcript.update((kind+str(box)+(repr(lam) if kind=='M' else '')).encode())
            row={'kind':kind,'box_unit_coordinates':[[str(a),str(b)] for a,b in box],'depth':depth}
            if kind=='M':row['lambda_exact']=str(Fraction(*lam.as_integer_ratio()))
            leaves.append(row)
            continue
        assert depth<100 and counts['checked']<100000
        middle=(l+h)*.5
        changes=[float(np.max(np.abs(np.diff(middle,axis=axis)))) for axis in range(3)]
        axis=int(np.argmax(changes))
        hl,hu,hr,hw=split_float_interval(l,h,axis)
        gl1,gh1,gl2,gh2=split_float_interval(gl,gh,axis)
        IH1,IH2=split_integer_upper(HI,axis)
        IG1,IG2=split_integer_upper(GI,axis)
        a,b=box[axis];m=(a+b)/2;box1=list(box);box2=list(box);box1[axis]=(a,m);box2[axis]=(m,b)
        stack.extend([(hr,hw,gl2,gh2,IH2,IG2,tuple(box2),depth+1),(hl,hu,gl1,gh1,IH1,IG1,tuple(box1),depth+1)])
        if counts['checked']%500==0:print('INTEGER_REPLAY',counts['checked'],'seconds',monotonic()-start,flush=True)
    assert covered==1
    for key,value in counts.items():assert value==expected[key],(key,value,expected[key])
    assert maxdepth==expected['max_depth']
    assert transcript.hexdigest()==expected['partition_transcript_sha256']
    leaves_raw=(json.dumps(leaves,indent=2)+'\n').encode();(HERE/'gram_replayed_leaves.json').write_bytes(leaves_raw)
    assert minmargin>Fraction(1,200000000)
    result={
        'status':'PASS: completed full-domain covering independently verified with integer upper arithmetic',
        'arithmetic':'Exact integer upper bounds on a fixed 2^-80 grid; exact integer dyadic multiplier acceptance; exact rational box geometry',
        'fixed_point_bits':BITS,
        'exact_cyclic_bernstein_sha256':audit['exact_cyclic_bernstein_sha256'],
        'all_4913_gram_initial_intervals_verified_exactly':True,
        **counts,
        'max_depth':maxdepth,
        'full_unit_cube_covered_exact':str(covered),
        'matches_original_partition_transcript':True,
        'partition_transcript_sha256':transcript.hexdigest(),
        'leaves_file':'gram_replayed_leaves.json',
        'leaves_sha256':hashlib.sha256(leaves_raw).hexdigest(),
        'rigorous_uniform_cyclic_negative_margin_exact':str(minmargin),
        'safe_strict_uniform_cyclic_upper_bound_exact':'-1/200000000',
        'margin_approximation_for_orientation_only':float(minmargin),
        'cyclic_gram_script_sha256':hashlib.sha256((CYCLIC/'bernstein_gram.py').read_bytes()).hexdigest(),
        'not_a_complete_count_or_density_theorem_by_itself':True,
        'runtime_seconds':monotonic()-start,
    }
    (HERE/'gram_integer_replay_result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
