"""Conservative exact cyclic decrement bounds for exterior graphs A, B, and E.

This helper supplies lower bounds on the *subtracted* correction. It does
not evaluate any base candidate or assert any global sign condition.
"""
from fractions import Fraction as Q

NEAR=(Q(1),Q(103,100))
MIDDLE=(Q(103,100),Q(3,2))
OUTER=(Q(3,2),Q(7,4))
ANTIPODAL_RADIAL=(Q(33,28),Q(7,4))
ANTIPODAL_ANGULAR=(Q(-1),Q(-7,8))
WEIGHTS=(Q(0),Q(0),Q(1,10000000))
CAPACITIES=(862,1102,1358)

def contained(box,interval,strict_left=False):
    lo,hi=box;a,b=interval
    return (a<lo if strict_left else a<=lo) and hi<=b

def quotient_interval(numerator,denominator):
    a,b=numerator;c,d=denominator
    assert a<=b and 0<c<=d
    vals=[a/c,a/d,b/c,b/d]
    return min(vals),max(vals)

def cosine_full(a,b,c):
    alo,ahi=a;blo,bhi=b;clo,chi=c
    assert 0<alo<=ahi and 0<blo<=bhi and 0<=clo<=chi
    return quotient_interval((alo*alo+blo*blo-chi*chi,ahi*ahi+bhi*bhi-clo*clo),(2*alo*blo,2*ahi*bhi))

def cosine_exterior(a,b,z):
    alo,ahi=a;blo,bhi=b;zlo,zhi=z
    assert Q(1)<=alo<=ahi and Q(1)<=blo<=bhi and 0<=zlo<=zhi<=1
    numlo=(alo+blo)**2-Q(49,16);numhi=(ahi+bhi)**2-Q(49,16)
    assert numlo>0
    lo,hi=quotient_interval((zlo*numlo,zhi*numhi),(2*alo*blo,2*ahi*bhi))
    return lo-1,hi-1

def mixed_band_indicator(a,b,t,far,angular,far_open_left=False):
    radial=(contained(a,NEAR) and contained(b,far,far_open_left)) or (contained(b,NEAR) and contained(a,far,far_open_left))
    return int(radial and contained(t,angular))

def cyclic_decrement_lower(kind,bounds,weights=WEIGHTS):
    """Rational lower bound valid throughout the given closed coordinate box.

    Full coordinates are triangle side lengths a,b,c. Exterior coordinates
    are a,b,z with t=-1+z*((a+b)^2-(7/4)^2)/(2ab), as in the retained verifier.
    A straddling support/angle interval contributes zero until subdivision.
    """
    assert len(bounds)==3 and all(Q(a)<=Q(b) for a,b in bounds)
    bounds=[tuple(map(Q,b)) for b in bounds]
    da,db,de=weights
    assert min(weights)>=0
    if kind=='outside':
        a,b,z=bounds;t=cosine_exterior(a,b,z)
        return da*mixed_band_indicator(a,b,t,MIDDLE,(Q(-4,5),Q(-1,2)),True)+db*mixed_band_indicator(a,b,t,OUTER,(Q(-19,20),Q(-3,4)))+de*int(contained(a,ANTIPODAL_RADIAL) and contained(b,ANTIPODAL_RADIAL) and contained(t,ANTIPODAL_ANGULAR))
    assert kind=='full'
    # All three graphs have selected distance strictly larger than 7/4.
    return Q(0)

def density_cost(weights=WEIGHTS):return sum((Q(c)*w for c,w in zip(CAPACITIES,weights)),Q(0))/8


def cyclic_decrement_upper(kind,bounds,weights=WEIGHTS):
    """Safe upper bound; at exact points the lower membership is exact."""
    assert kind in ('full','outside') and len(bounds)==3 and min(weights)>=0
    if kind=='full':return Q(0)
    if all(Q(lo)==Q(hi) for lo,hi in bounds):
        return cyclic_decrement_lower(kind,bounds,weights)
    return sum(weights,Q(0))
