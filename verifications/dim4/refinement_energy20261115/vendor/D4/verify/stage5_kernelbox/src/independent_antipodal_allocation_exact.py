"""Exact local data and box bounds for an independent antipodal allocation.

This separate helper covers the inner_B4 profile, support [1,7/4], and
zero displacement of the reflected test point.  It does not certify a
global kernel.  In particular, the caller must include weight*u in the
complete unlifted pair residual and freshly prove that residual <= 0
before constructing its count-gauge budget.

The ordinary capacity follows from the finite-set variance identity.
The optional classical improvement has the elementary capacity proof in
copositive/CLASSICAL_SHIFTED_ANTIPODAL_CAPACITY.md, Section 1.  This helper
uses only its unshifted formula; movable-center records are rejected.
All returned endpoints are rational, with integer-square-root enclosures
for the classical curved piece.  No floating point is used.
"""
from fractions import Fraction as Q
from itertools import product
from math import comb, isqrt


PROFILES = ('ordinary', 'classical_improvement')
RADIUS = Q(7,4)
INNER_INTERVAL = (Q(1), Q(3,2))
OUTER_INTERVAL = (Q(3,2), Q(7,4))


def rational_sqrt_bounds(value, bits=96):
    """Outward dyadic bounds, exact even at a rational square endpoint."""
    value = Q(value)
    assert value >= 0 and isinstance(bits,int) and bits >= 1
    denominator = 1 << bits
    integer = isqrt(value.numerator*denominator*denominator//value.denominator)
    lower = Q(integer,denominator)
    if integer*integer*value.denominator == value.numerator*denominator*denominator:
        return lower,lower
    return lower,Q(integer+1,denominator)


def capacity_point_bounds(squared_distance, capacity_profile='ordinary', bits=96):
    """Rational enclosure of the exact capacity at one rational square norm."""
    s = Q(squared_distance)
    assert s >= 0 and capacity_profile in PROFILES
    if s >= Q(1,2):
        return Q(0),Q(0)
    if capacity_profile == 'ordinary' or s >= Q(1,4):
        value = 1-2*s
        return value,value
    lower_root,upper_root = rational_sqrt_bounds(s,bits)
    # s<1/4 and the dyadic endpoint 1/2 imply upper_root<=1/2.
    assert 0 <= lower_root <= upper_root <= Q(1,2)
    lower = min(Q(1),2*(1-upper_root)**2)
    upper = min(Q(1),2*(1-lower_root)**2)
    assert 0 <= lower <= upper <= 1
    return lower,upper


def squared_distance_bounds(radius, pivot, coordinate, cutoff=RADIUS):
    """Sharp rational range of (r-s)^2+z*((r+s)^2-R^2) on a box.

    The expression increases in z.  At fixed z it is convex in (r,s),
    so its maximum uses corners.  Its minimum uses corners and the
    stationary point on each of four edges.  For z>0 there is no interior
    stationary point with r,s>=1; for z=0 its zero ridge intersects an
    edge whenever the two radius intervals overlap.
    """
    rlo,rhi = map(Q,radius)
    slo,shi = map(Q,pivot)
    zlo,zhi = map(Q,coordinate)
    cutoff = Q(cutoff)
    assert 1 <= rlo <= rhi <= cutoff <= 2
    assert 1 <= slo <= shi <= cutoff and 0 <= zlo <= zhi <= 1
    def value(r,s,z):
        return (r-s)**2+z*((r+s)**2-cutoff*cutoff)
    lower_candidates = [value(r,s,zlo) for r,s in product((rlo,rhi),(slo,shi))]
    ratio = (1-zlo)/(1+zlo)
    for s in (slo,shi):
        r = ratio*s
        if rlo <= r <= rhi:
            lower_candidates.append(value(r,s,zlo))
    for r in (rlo,rhi):
        s = ratio*r
        if slo <= s <= shi:
            lower_candidates.append(value(r,s,zlo))
    lower = min(lower_candidates)
    upper = max(value(r,s,zhi) for r,s in product((rlo,rhi),(slo,shi)))
    assert 0 <= lower <= upper
    return lower,upper


def capacity_box_bounds(radius, pivot, coordinate, capacity_profile='ordinary',
                        cutoff=RADIUS, bits=96):
    """Use the monotone capacity and the exact squared-distance box range."""
    lower_squared,upper_squared = squared_distance_bounds(radius,pivot,coordinate,cutoff)
    lower = capacity_point_bounds(upper_squared,capacity_profile,bits)[0]
    upper = capacity_point_bounds(lower_squared,capacity_profile,bits)[1]
    assert 0 <= lower <= upper <= 1
    return lower,upper


def embed_inner_b4(degree,basis_intervals):
    """Exact degree elevation into broad-interval regional Bernstein blocks."""
    assert isinstance(degree,int) and degree >= 4
    result = []
    for interval in basis_intervals:
        interval = tuple(map(Q,interval))
        assert interval in (INNER_INTERVAL,OUTER_INTERVAL)
        if interval == INNER_INTERVAL:
            result.extend(Q(comb(index,4),comb(degree,4)) if index >= 4 else Q(0)
                          for index in range(degree+1))
        else:
            result.extend([Q(0)]*(degree+1))
    return result


def validate_allocations(records,core,basis_intervals):
    """Validate this narrow exact allocation schema, including physical w.

    The caller has already parsed the primary kernel's rational partition.
    Returned weights multiply local energies directly, with no hidden
    capacity normalization.  Their origin constant is zero.
    """
    assert core.NC == core.DEG+1 and core.DEG >= 4
    assert core.NS == 7 and len(core.BREAKS) == 8
    assert core.BREAKS[0] == 1 and tuple(map(Q,core.BREAKS[-2:])) == (Q(3,2),RADIUS)
    assert len(basis_intervals) == core.NS
    for region,(a,b) in enumerate(basis_intervals):
        a,b = Q(a),Q(b)
        assert a <= core.BREAKS[region] < core.BREAKS[region+1] <= b
        assert (a,b) == (INNER_INTERVAL if region < 6 else OUTER_INTERVAL)
    expected_w = embed_inner_b4(core.DEG,basis_intervals)
    assert len(expected_w) == core.NC*core.NS
    expected_profile = tuple([Q(0)]*4+[Q(1)]+[Q(0)]*10)
    answer = []
    for source in records:
        assert source['kind'] == 'independent_antipodal_capacity'
        assert source['profile'] == 'inner_B4' and source['radial_degree'] == 4
        assert source['capacity_profile'] in PROFILES
        assert tuple(map(Q,source['support'])) == (Q(1),RADIUS)
        assert tuple(tuple(map(Q,row)) for row in source['broad_intervals']) == (
            INNER_INTERVAL,OUTER_INTERVAL)
        assert tuple(map(Q,source['profile_coefficients'])) == expected_profile
        assert tuple(map(Q,source['indicator_radii'])) == tuple(map(Q,core.BREAKS[1:6]))
        assert source['weight_normalization'] == 'actual_local_energy_multiplier'
        assert Q(source['origin_constant']) == 0
        assert Q(source.get('shift','0')) == 0 and Q(source.get('delta','0')) == 0
        assert 'pivot_lower_cutoff' not in source
        weight = Q(source['weight'])
        w = list(map(Q,source['w']))
        assert weight >= 0 and w == expected_w
        answer.append(dict(source,weight=weight,w=w,
            capacity_validation='FINITE_UNIT_SEPARATED_TRANSLATED_CAPACITY',
            capacity_proof=('FINITE_VARIANCE_IDENTITY' if source['capacity_profile']=='ordinary'
                            else 'CLASSICAL_BLICHFELDT_IMPROVEMENT_ELEMENTARY_PROOF')))
    return answer


def profile_bounds(radius,region,basis_intervals):
    """Bounds for u, preserving the actual-left and right-closure values."""
    lo,hi = map(Q,radius)
    interval = tuple(map(Q,basis_intervals[region]))
    assert interval[0] <= lo <= hi <= interval[1]
    if interval == OUTER_INTERVAL:
        return Q(0),Q(0)
    assert interval == INNER_INTERVAL
    return (2*lo-2)**4,(2*hi-2)**4


def correction_bounds(kind,bounds,regions,allocations,core,basis_intervals,bits=96):
    """Bound the nonnegative ordinary capacity quantity subtracted from Gamma.

    This is only the independent antipodal correction.  The complete pair
    origin shift and the resulting count budget belong to the caller's
    existing exact pair and count-gauge constructions.
    """
    assert kind in ('full','outside') and len(bounds) == 3
    assert len(regions) == (3 if kind=='full' else 2)
    assert Q(core.BREAKS[-1]) == RADIUS
    for axis,region in enumerate(regions):
        lo,hi = map(Q,bounds[axis])
        assert core.BREAKS[region] <= lo <= hi <= core.BREAKS[region+1]
    if kind == 'full':
        assert 4-RADIUS*RADIUS > Q(1,2)
        return Q(0),Q(0)
    lower = upper = Q(0)
    u0 = profile_bounds(bounds[0],regions[0],basis_intervals)
    u1 = profile_bounds(bounds[1],regions[1],basis_intervals)
    for allocation in allocations:
        assert allocation['capacity_validation'] == 'FINITE_UNIT_SEPARATED_TRANSLATED_CAPACITY'
        weight = Q(allocation['weight'])
        assert weight >= 0
        phi = capacity_box_bounds(*bounds,capacity_profile=allocation['capacity_profile'],
                                  cutoff=RADIUS,bits=bits)
        lower += weight*(u0[0]+u1[0])*phi[0]/2
        upper += weight*(u0[1]+u1[1])*phi[1]/2
    assert 0 <= lower <= upper
    return lower,upper
