"""Separate exact signed-exponent int64 handoff, preserving the old API.

For interval integers L,U at input scale S, the physical coefficient lies
in [L/S,U/S]. For any integer q, including q<0, outward rounding at scale
sigma=2**q gives int64 endpoints l,u with l/sigma <= coefficient <= u/sigma.
The resulting positive Fraction scale MUST also be used in every polynomial
and correction comparison. This function makes no sign decision.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib

import numpy as np
import fixed_binary_bernstein_v3 as base

BASE_SHA256='dee3a3fafa6ba7cb811bae7d2568b8633ed68eed9702155132cac79f70075f07'
if hashlib.sha256(Path(base.__file__).read_bytes()).hexdigest()!=BASE_SHA256:
    raise ValueError('Interval backend dependency hash mismatch.')


def to_int64_signed_exponent(low, high, scale, *, maximum_exponent=50,
                             minimum_exponent=-4096):
    low,high=base._tensor(low,'lower'),base._tensor(high,'upper')
    if low.shape!=high.shape or any(int(l)>int(u) for l,u in zip(low.flat,high.flat)):
        raise ValueError('Mismatched or inverted interval endpoints.')
    scale=base._integer(scale,'scale')
    if scale<=0 or scale&(scale-1):
        raise ValueError('Input scale must be a positive power of two.')
    qmax=base._integer(maximum_exponent,'maximum_exponent')
    qmin=base._integer(minimum_exponent,'minimum_exponent')
    if not -4096<=qmin<=qmax<=50:
        raise ValueError('Require -4096 <= minimum <= maximum <= 50.')
    minimum,maximum=min(map(int,low.flat)),max(map(int,high.flat))
    limit=1<<59
    for q in range(qmax,qmin-1,-1):
        numerator,denominator=(1<<q,scale) if q>=0 else (1,scale*(1<<(-q)))
        lo_min=minimum*numerator//denominator
        hi_max=base._ceil_div(maximum*numerator,denominator)
        if max(abs(lo_min),abs(hi_max))<limit:
            lo=base._round_ratio(low,numerator,denominator,False)
            hi=base._round_ratio(high,numerator,denominator,True)
            if not all(-limit<int(v)<limit for a in (lo,hi) for v in a.flat):
                raise OverflowError('Signed-exponent int64 overflow guard failed.')
            output_scale=Q(1<<q) if q>=0 else Q(1,1<<(-q))
            return lo.astype(np.int64),hi.astype(np.int64),output_scale,q
    raise OverflowError('Dynamic range exceeds the explicit signed-exponent limit.')
