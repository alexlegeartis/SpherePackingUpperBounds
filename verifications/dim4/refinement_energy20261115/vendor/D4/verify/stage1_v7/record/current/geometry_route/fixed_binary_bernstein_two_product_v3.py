"""Two-product fixed-binary enclosure; separate from the four-product base.

If T_L=floor(S*T), then 0<=S*T-T_L<=1 entrywise. With L<=S*C<=U,
put P=T_L*L and W=T_L*(U-L). For every output column, lower numerator
P+sum(min(L,0)) and upper numerator P+W+sum(max(U,0)) enclose S^2*T*C.
Column sums are broadcast to all output rows. Outward division by S keeps
the original fixed scale. No coefficient/sign conclusion uses sampling.
"""
from pathlib import Path
import hashlib
import time

import numpy as np
import fixed_binary_bernstein_v3 as base

BASE_SHA256='dee3a3fafa6ba7cb811bae7d2568b8633ed68eed9702155132cac79f70075f07'
if hashlib.sha256(Path(base.__file__).read_bytes()).hexdigest()!=BASE_SHA256:
    raise ValueError('Four-product base dependency hash mismatch.')


def _product(TL,L,U,S,column_block):
    # Width is nonnegative, so its product often has much smaller integer
    # entries than the signed central product. Both are exact over Z.
    P=base._reference.integer_matmul(TL,L,column_block)
    W=base._reference.integer_matmul(TL,U-L,column_block)
    negative=np.minimum(L,0).sum(axis=0)
    positive=np.maximum(U,0).sum(axis=0)
    lower=base._round_ratio(P+negative,1,S,False)
    upper=base._round_ratio(P+W+positive,1,S,True)
    if any(int(a)>int(b) for a,b in zip(lower.flat,upper.flat)):
        raise ArithmeticError('Inverted two-product interval.')
    return lower,upper


def to_bernstein_interval(C,denominator,bounds,*,precision=384,
                          column_block=1024,progress=None):
    precision=base._integer(precision,'precision')
    if not 1<=precision<=4096:
        raise ValueError('precision must lie in [1,4096].')
    column_block=base._integer(column_block,'column_block')
    if not 1<=column_block<=65536:
        raise ValueError('column_block must lie in [1,65536].')
    denominator=base._integer(denominator,'denominator')
    if denominator<=0:
        raise ValueError('denominator must be positive.')
    C=base._tensor(C,'C')
    if not isinstance(bounds,(list,tuple)) or len(bounds)!=C.ndim:
        raise ValueError('One exact interval per tensor axis is required.')
    parsed=[]
    for b in bounds:
        if not isinstance(b,(list,tuple)) or len(b)!=2:
            raise ValueError('Each bound must be a pair.')
        left,right=(base._exact_rational(v,'bound') for v in b)
        if not 0<=left<right:
            raise ValueError('Every axis requires 0 <= left < right.')
        parsed.append((left,right))
    if progress is not None and not callable(progress):
        raise TypeError('progress must be callable or None.')
    S=1<<precision
    low=base._round_ratio(C,S,denominator,False)
    high=base._round_ratio(C,S,denominator,True)
    del C
    for axis,(left,right) in enumerate(parsed):
        started=time.monotonic();n=low.shape[axis]-1
        TL,_=base._rounded_transform(n,left,right-left,precision)
        L,U=np.moveaxis(low,axis,0),np.moveaxis(high,axis,0)
        shape=L.shape
        lo,hi=_product(TL,L.reshape(n+1,-1),U.reshape(n+1,-1),S,column_block)
        low=np.moveaxis(lo.reshape(shape),0,axis)
        high=np.moveaxis(hi.reshape(shape),0,axis)
        if progress is not None:
            progress(dict(axis=axis,degree=n,seconds=time.monotonic()-started,precision=precision,
                maximum_endpoint_bits=max(abs(int(v)).bit_length() for a in (low,high) for v in a.flat),
                maximum_width_units=max(int(b)-int(a) for a,b in zip(low.flat,high.flat))))
    return low,high,S
