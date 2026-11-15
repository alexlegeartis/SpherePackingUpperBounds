"""Exact rational certificates for negative-pole pointed tent corrections.

This is a separate extension; accepted pointed-ball checker files are not
modified.  A tent record has rational fields d, capacity, weight, w and
optional sigma=-1.  Its unnormalized coefficient is beta=weight:

    h(r)=capacity*beta*w(r)/2,
    K(x,y)=-beta*(w(s)*phi(x+d*yhat)+w(r)*phi(y+d*xhat))/2,
    phi(z)=max(0,1-2*|z|^2).

Both nonorigin arguments are restricted to the core annulus.  In particular
K(x,x)=0 and the pair contribution is capacity*beta*w(r).  Kernel origin
mass and its other local-energy summands are verified by the calling code.
This helper only validates the tent summands and the cyclic inequalities.
"""
from fractions import Fraction as Q
import itertools
import time
import numpy as np
from interval_certificate import radial_bounds, sqrt_bounds
from interval_certificate_full_closed import correction_bounds as cap_correction_bounds


def validate_capacity(distance, capacity):
    """Prove a rational nonorigin tent-capacity bound without real arithmetic."""
    d, c = Q(distance), Q(capacity)
    assert d > 0 and c > 0
    if c >= 1:
        return 'UNIVERSAL_PAIR_VARIANCE_TENT_BOUND'
    assert d*d <= Q(1, 2)
    if c >= 2*d*d:
        return 'UNIVERSAL_TENT_BOUND_MINUS_ORIGIN_TENT'
    # At most four nonorigin points have positive tent weight: together
    # with the origin they form a strictly obtuse set in R4.  Every possible
    # cardinality n has capacity F_n(d) below.  Since the left hand side of
    # the unsquared comparison is positive, these squared tests are exact.
    for n in range(1, 5):
        lhs = c+n*(1+2*d*d)
        rhs_squared = 8*d*d*n*(n+1)
        assert lhs > 0 and lhs*lhs >= rhs_squared, (n, d, c)
    return 'R4_FOUR_SIMPLEX_TENT_CAPACITY_ENVELOPE'


def validate_tents(records, core, basis_intervals):
    """Parse nonnegative physical Bernstein rays and validate their capacities."""
    assert core.NC == core.DEG+1 and core.DEG >= 1
    assert len(basis_intervals) == core.NS
    assert core.BREAKS[0] == 1 and core.BREAKS[-1] <= 2
    for i, (a, b) in enumerate(basis_intervals):
        assert a <= core.BREAKS[i] < core.BREAKS[i+1] <= b
    tents = []
    for record in records:
        tent = {key: Q(record[key]) for key in ['d', 'capacity', 'weight']}
        tent['sigma'] = Q(record.get('sigma', -1))
        tent['w'] = list(map(Q, record['w']))
        assert tent['sigma'] == -1 and tent['weight'] >= 0
        assert len(tent['w']) == core.NC*core.NS
        assert all(coefficient >= 0 for coefficient in tent['w'])
        # No truncation convention is implicit in a tent record.  Radial
        # step functions must be represented in the actual branch weights.
        assert 'pivot_lower_cutoff' not in record
        tent['capacity_lemma'] = validate_capacity(tent['d'], tent['capacity'])
        tents.append(tent)
    return tents


def tent_bounds(radius, cosine, tent):
    """Exact range of the polynomial hinge on a radius/cosine rectangle."""
    lo, hi = radius
    tlo, thi = cosine
    d = tent['d']
    assert d > 0 and 0 < lo <= hi and -1 <= tlo <= thi <= 1
    value = lambda r,t: 1-2*(r*r+d*d)-4*d*r*t
    # The expression decreases in t and is concave in r.
    lower = min(value(lo,thi), value(hi,thi))
    critical = min(hi,max(lo,-d*tlo))
    upper = value(critical,tlo)
    return max(Q(0),lower), max(Q(0),upper)


def full_tent_bounds(radius, pivot, third, tent):
    """Sharp box enclosure using the triangle formula, before taking the hinge.

P=max(0,1-2(r^2+d^2)-2d*(r^2+s^2-z^2)/s).  It decreases in r,
increases in z, and at fixed (r,z) is monotone in s or has one interior
maximum.  Its lower bound uses corners only.  An interior upper extremum
contains sqrt(r^2-z^2), which receives an outward rational enclosure.
"""
    rlo,rhi = radius
    slo,shi = pivot
    zlo,zhi = third
    d = tent['d']
    assert min(rlo,slo,zlo)>0 and rlo<=rhi and slo<=shi and zlo<=zhi
    value = lambda r,s,z: 1-2*(r*r+d*d)-2*d*(r*r+s*s-z*z)/s
    lower = min(value(rhi,s,zlo) for s in pivot)
    upper = max(value(rlo,s,zhi) for s in pivot)
    square = rlo*rlo-zhi*zhi
    if square > 0 and slo*slo <= square <= shi*shi:
        upper = max(upper,1-2*(rlo*rlo+d*d)-4*d*sqrt_bounds(square)[0])
    return max(Q(0),lower), max(Q(0),upper)


def outside_tent_bounds(radius, pivot, coordinate, cutoff, tent):
    """Exact box range in the exterior branch's normalized third coordinate.

t=-1+z*((r+s)^2-R^2)/(2rs), 0<=z<=1.  The tent decreases in
z and in s on 1<=r,s<=R<=2; it is concave in r.  Both extrema below
therefore follow from rational endpoint and stationary-point evaluations.
"""
    rlo,rhi = radius
    slo,shi = pivot
    zlo,zhi = coordinate
    d = tent['d']
    assert 1 <= rlo <= rhi <= cutoff <= 2
    assert 1 <= slo <= shi <= cutoff and 0 <= zlo <= zhi <= 1
    value = lambda r,s,z: 1-2*(r-d)**2-2*d*z*((r+s)**2-cutoff*cutoff)/s
    lower = min(value(r,shi,zhi) for r in radius)
    critical = d*(1-zlo)/(1+d*zlo/slo)
    critical = min(rhi,max(rlo,critical))
    upper = value(critical,slo,zlo)
    return max(Q(0),lower), max(Q(0),upper)


def full_tents_vanish(cutoff):
    """Exact test for R^2<=2+sqrt(2), where all negative-pole tents vanish.

Positive tent weight implies x.yhat< -sqrt(|x|^2-1/2), so the third
packing distance has square strictly greater than 2+sqrt(2).  The strict
support handles the equality endpoint without any convention about balls.
"""
    excess = Q(cutoff)**2-2
    return excess<=0 or excess*excess<=2


def correction_bounds(kind, bounds, regions, tents, core, basis_intervals):
    """Enclose the positive quantity subtracted from the cyclic expression."""
    assert kind in ['full','outside'] and len(bounds)==3
    assert len(regions)==(3 if kind=='full' else 2)
    if kind=='full' and full_tents_vanish(core.BREAKS[-1]):
        assert all(1<=lo<=hi<=core.BREAKS[-1] for lo,hi in bounds)
        return Q(0),Q(0)
    lower = upper = Q(0)
    for tent in tents:
        w = [radial_bounds(bounds[i],regions[i],tent,basis_intervals,core.DEG)
             for i in range(len(regions))]
        pairs = [(0,1,2),(0,2,1),(1,2,0)] if kind=='full' else [(0,1,2)]
        for i,j,k in pairs:
            if kind=='full':
                pi = full_tent_bounds(bounds[i],bounds[j],bounds[k],tent)
                pj = full_tent_bounds(bounds[j],bounds[i],bounds[k],tent)
            else:
                pi = outside_tent_bounds(bounds[i],bounds[j],bounds[k],core.BREAKS[-1],tent)
                pj = outside_tent_bounds(bounds[j],bounds[i],bounds[k],core.BREAKS[-1],tent)
            beta = tent['weight']/2
            lower += beta*(w[j][0]*pi[0]+w[i][0]*pj[0])
            upper += beta*(w[j][1]*pi[1]+w[i][1]*pj[1])
    return lower,upper


def check(core, kind, regions, eta, tents, basis_intervals,
          node_limit=180000, depth_limit=66, caps=()):
    """Prove a cyclic branch including tent and optional pointed-ball rays.

The caller must first call validate_tents, and independently validate any
pointed caps using interval_certificate_full.validate_caps.  All additional
terms have nonnegative radial coefficients and nonpositive off-diagonal
kernel entries.  This routine supplies no assertion about pair constraints,
origin mass, the scalar function, or unvisited branches.
"""
    started = time.time()
    assert kind in ['full','outside']
    C,denominator,domain = core.build(kind,regions,eta=eta)
    C,denominator = core.to_bernstein(C,denominator,domain)
    biggest = max(abs(int(v)) for v in C.flat)
    precision = min(50,58-(biggest.bit_length()-denominator.bit_length()+1))
    if precision < 0:
        return dict(status='UNKNOWN',reason='coefficient dynamic range')
    scale = 1 << precision
    low = np.array([int(v)*scale//denominator for v in C.flat],dtype=np.int64).reshape(C.shape)
    high = np.array([-((-int(v)*scale)//denominator) for v in C.flat],dtype=np.int64).reshape(C.shape)
    assert max(abs(low).max(),abs(high).max()) < 2**59
    del C
    angular_degree = core.gram()[2]

    def factor(point):
        value = 2*point[0]*point[1]
        if kind=='full':
            value *= point[2]
        return value**angular_degree

    def all_corrections(bounds):
        tent_lo,tent_hi = correction_bounds(kind,bounds,regions,tents,core,basis_intervals)
        cap_lo,cap_hi = (cap_correction_bounds(kind,bounds,regions,caps,core,basis_intervals)
                         if caps else (Q(0),Q(0)))
        return tent_lo+cap_lo,tent_hi+cap_hi,tent_lo,cap_lo

    stack = [(low,high,((0,1),)*3,(0,)*3)]
    nodes = leaves = correction_leaves = tent_dependent_leaves = maximum_depth = 0
    failures = []
    unknown = False
    while stack:
        lo,hi,box,levels = stack.pop()
        nodes += 1
        maximum_depth = max(maximum_depth,sum(levels))
        maximum = int(hi.max())
        if maximum <= 0:
            leaves += 1
            continue
        bounds = [(left+(right-left)*Q(x,2**level),left+(right-left)*Q(y,2**level))
                  for (left,right),(x,y),level in zip(domain,box,levels)]
        correction_lower,_,_,cap_lower = all_corrections(bounds)
        factor_lower = factor([b[0] for b in bounds])
        if Q(maximum,scale) <= factor_lower*correction_lower:
            leaves += 1
            correction_leaves += 1
            tent_dependent_leaves += int(Q(maximum,scale)>factor_lower*cap_lower)
            continue
        failure = None
        for corner in itertools.product((0,1),repeat=3):
            index = tuple(-1 if side else 0 for side in corner)
            base_lower = int(lo[index])
            if base_lower <= 0:
                continue
            point = [bounds[i][side] for i,side in enumerate(corner)]
            _,correction_upper,_,_ = all_corrections([(p,p) for p in point])
            if Q(base_lower,scale)>factor(point)*correction_upper:
                failure = dict(point=list(map(str,point)),base_integer_lower=base_lower,
                               correction_upper=str(correction_upper))
                break
        if failure:
            failures.append(failure)
            break
        if nodes>=node_limit or sum(levels)>=depth_limit:
            unknown = True
            break
        axis = min(range(3),key=lambda j:levels[j])
        ll,lr = core.split(lo,axis,False)
        ul,ur = core.split(hi,axis,True)
        nextlevels = list(levels)
        nextlevels[axis] += 1
        leftbox,rightbox = list(box),list(box)
        a,b = box[axis]
        leftbox[axis],rightbox[axis] = (2*a,a+b),(a+b,2*b)
        stack.extend([(lr,ur,tuple(rightbox),tuple(nextlevels)),
                      (ll,ul,tuple(leftbox),tuple(nextlevels))])
    status = 'VIOLATION' if failures else ('UNKNOWN' if unknown or stack else 'CERTIFIED')
    record = dict(status=status,kind=kind,regions=list(regions),nodes=nodes,leaves=leaves,
                  correction_leaves=correction_leaves,tent_dependent_leaves=tent_dependent_leaves,
                  maxdepth=maximum_depth,binary_scale=precision,failures=failures,
                  seconds=time.time()-started,
                  pointed_helper='EXACT_NEGATIVE_POLE_TENT_WITH_OPTIONAL_POINTED_CAPS')
    print('EXACT POINTED TENT',record,flush=True)
    return record
