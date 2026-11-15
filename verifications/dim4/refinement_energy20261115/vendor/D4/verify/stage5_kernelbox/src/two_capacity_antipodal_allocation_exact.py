"""Exact antipodal allocation bounds including a rational two-level capacity.

The frozen single-capacity helper supplies the already audited physical
profile, box-distance geometry and classical capacity. This separate
module adds only the exact step function

  1 if sqrt(s)<a; 1/2 if a<=sqrt(s)<=b; 0 if sqrt(s)>b,

where a=8453/20000, b=11547/20000, a+b=1 and 3*b*b<1.
A closed b-ball contains at most two unit-separated points, by variance.
A point strictly inside a excludes all other support points. Thus every
translate has capacity at most one, including both prescribed endpoints.

This is a local helper, not a global packing certificate. All independent
allocation pair increments must enter the caller's complete unlifted q
before its separate count-gauge budget is formed and verified.
"""
from fractions import Fraction as Q
import independent_antipodal_allocation_exact as frozen


RADIUS = frozen.RADIUS
INNER_INTERVAL = frozen.INNER_INTERVAL
OUTER_INTERVAL = frozen.OUTER_INTERVAL
INNER = Q(8453,20000)
OUTER = Q(11547,20000)
assert INNER > 0 and INNER+OUTER == 1 and 3*OUTER*OUTER < 1
assert OUTER*OUTER < Q(1,2)
PROFILES = frozen.PROFILES+('two_level',)

# Exact, unchanged functionality of the frozen single-capacity helper.
rational_sqrt_bounds = frozen.rational_sqrt_bounds
squared_distance_bounds = frozen.squared_distance_bounds
embed_inner_b4 = frozen.embed_inner_b4
profile_bounds = frozen.profile_bounds


def capacity_point_bounds(squared_distance,capacity_profile='ordinary',bits=96):
    """Exact endpoint conventions for the new step; old profiles delegate."""
    s = Q(squared_distance)
    assert s >= 0 and capacity_profile in PROFILES
    if capacity_profile != 'two_level':
        return frozen.capacity_point_bounds(s,capacity_profile,bits)
    if s < INNER*INNER:
        value = Q(1)
    elif s <= OUTER*OUTER:
        value = Q(1,2)
    else:
        value = Q(0)
    return value,value


def capacity_box_bounds(radius,pivot,coordinate,capacity_profile='ordinary',
                        cutoff=RADIUS,bits=96):
    """Monotone capacity applied to exact squared-distance box endpoints."""
    lower_squared,upper_squared = squared_distance_bounds(radius,pivot,coordinate,cutoff)
    lower = capacity_point_bounds(upper_squared,capacity_profile,bits)[0]
    upper = capacity_point_bounds(lower_squared,capacity_profile,bits)[1]
    assert 0 <= lower <= upper <= 1
    return lower,upper


def validate_allocations(records,core,basis_intervals):
    """Reuse every old schema gate and add only the fixed rational capacity.

    For the new record, temporarily select the ordinary profile solely to
    apply the frozen radial, support, weight and zero-origin validations.
    The returned record restores the step profile and records its distinct
    independently justified capacity proof. No numerical formula is read.
    """
    assert isinstance(records,list)
    answer = []
    for source in records:
        assert source['capacity_profile'] in PROFILES
        if source['capacity_profile'] != 'two_level':
            assert not any(key in source for key in
                           ('inner_radius','outer_radius','inner_boundary','outer_boundary'))
            answer.extend(frozen.validate_allocations([source],core,basis_intervals))
            continue
        assert Q(source['inner_radius']) == INNER and Q(source['outer_radius']) == OUTER
        assert source['inner_boundary'] == 'strict' and source['outer_boundary'] == 'closed'
        assert Q(source.get('inner_weight','1')) == 1
        assert Q(source.get('outer_weight','1/2')) == Q(1,2)
        delegated = dict(source,capacity_profile='ordinary')
        checked = frozen.validate_allocations([delegated],core,basis_intervals)[0]
        checked.update(capacity_profile='two_level',
            capacity_proof='CLOSED_B_BALL_HAS_AT_MOST_TWO_POINTS_AND_STRICT_A_POINT_EXCLUDES_ALL_OTHERS',
            exact_inner_radius=str(INNER),exact_outer_radius=str(OUTER),
            exact_three_point_variance_comparison=str(3*OUTER*OUTER))
        answer.append(checked)
    return answer


def correction_bounds(kind,bounds,regions,allocations,core,basis_intervals,bits=96):
    """Bounds on the nonnegative ordinary quantity subtracted from Gamma."""
    assert kind in ('full','outside') and len(bounds) == 3
    assert len(regions) == (3 if kind=='full' else 2)
    assert Q(core.BREAKS[-1]) == RADIUS
    for axis,region in enumerate(regions):
        lo,hi = map(Q,bounds[axis])
        assert core.BREAKS[region] <= lo <= hi <= core.BREAKS[region+1]
    if kind == 'full':
        # All old profiles vanish at squared norm >=1/2; the new one
        # vanishes already above b^2<1/2. This includes all closures.
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
