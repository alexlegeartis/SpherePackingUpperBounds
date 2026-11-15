"""Exact pointed correction bounds on full and exterior triangle branches.

This extends the separately audited general probability helper without
changing either accepted v1 checker.  A numerical optimizer is never used
here.  The full triangle cosine extrema are enclosed over the entire box;
the subtracted correction is compared against an outward integer Bernstein
bound for the positive homogenization of the polynomial kernel.
"""
from fractions import Fraction as Q
import itertools
import time
import numpy as np
from interval_certificate import sqrt_bounds, radial_bounds
from interval_certificate_general_closed import probability_bounds
from interval_certificate_general_closed import validate_caps as validate_general_caps
from interval_certificate_general_closed import correction_bounds as exterior_correction_bounds


def validate_caps(records, core, basis_intervals):
    """Validate ordinary probes and the exact onset-truncation contract."""
    caps = validate_general_caps(records, core, basis_intervals)
    for record, cap in zip(records, caps):
        if 'pivot_lower_cutoff' not in record:
            continue
        cutoff = Q(record['pivot_lower_cutoff'])
        assert cap['sigma'] == 1 and cap['M'] >= 1
        assert cap['d']**2 <= cap['a2'] < Q(1, 2)
        assert 1 <= cutoff < core.BREAKS[-1]
        onset = (1+cap['d']**2-cap['a2'])/cap['d']
        assert cutoff < onset, 'Closed lower piece must lie strictly below the geometric onset.'
        cap['pivot_lower_cutoff'] = cutoff
        cap['pivot_truncation_proof'] = 'DISTINCT_POINT_DISTANCE_SQUARED_AT_LEAST_1_PLUS_D_SQUARED_MINUS_D_TIMES_PIVOT_RADIUS'
    return caps


def angle_bounds(a, b, c):
    """Enclose (a^2+b^2-c^2)/(2ab) on a triangle side-length box.

    Dependence on c is decreasing.  At fixed c>0 there is no stationary
    point in the interior of the (a,b) rectangle: its two derivative
    equations would imply c=0.  On an edge there is at most one critical
    point, a minimum at sqrt(fixed_side^2-c^2).  Thus maxima occur at
    corners; minima occur at corners or the listed edge critical points.
    Every candidate square root receives an outward rational enclosure.
    """
    assert 0 < a[0] <= a[1] and 0 < b[0] <= b[1] and 0 < c[0] <= c[1]
    assert a[0]+b[0] >= c[1] and a[0]+c[0] >= b[1] and b[0]+c[0] >= a[1]
    cosine = lambda r, s, z: (r*r+s*s-z*z)/(2*r*s)
    lower = min(cosine(r, s, c[1]) for r, s in itertools.product(a, b))
    upper = max(cosine(r, s, c[0]) for r, s in itertools.product(a, b))
    for fixed, variable in [(a, b), (b, a)]:
        for side in fixed:
            square = side*side-c[1]*c[1]
            if square > 0 and variable[0]**2 <= square <= variable[1]**2:
                lower = min(lower, sqrt_bounds(square)[0]/side)
    lower, upper = max(Q(-1), lower), min(Q(1), upper)
    assert -1 <= lower <= upper <= 1
    return lower, upper


def correction_bounds(kind, bounds, regions, caps, core, basis_intervals):
    # A validated pivot truncation does not change any distinct-point kernel
    # value on a feasible packing triangle.  Hence untruncated radial weights
    # give the exact same cyclic correction, without a new radial partition.
    if kind == 'outside':
        return exterior_correction_bounds(bounds, regions, caps, core, basis_intervals)
    assert kind == 'full' and len(bounds) == len(regions) == 3
    assert all(1 <= lo <= hi <= core.BREAKS[-1] <= 2 for lo, hi in bounds)
    angles = {(i, j): angle_bounds(bounds[i], bounds[j], bounds[k])
              for i, j, k in [(0, 1, 2), (0, 2, 1), (1, 2, 0)]}
    lower = upper = Q(0)
    for cap in caps:
        w = [radial_bounds(bounds[i], regions[i], cap, basis_intervals, core.DEG)
             for i in range(3)]
        for i, j in angles:
            pi = probability_bounds(bounds[i], angles[i, j], cap)
            pj = probability_bounds(bounds[j], angles[i, j], cap)
            lower += cap['weight']*(w[j][0]*pi[0]+w[i][0]*pj[0])/2
            upper += cap['weight']*(w[j][1]*pi[1]+w[i][1]*pj[1])/2
    return lower, upper


def check(core, kind, regions, eta, caps, basis_intervals,
          node_limit=180000, depth_limit=66):
    """Prove the full cyclic expression including all pointed corrections."""
    started = time.time()
    assert kind in ['full', 'outside']
    C, denominator, domain = core.build(kind, regions, eta=eta)
    C, denominator = core.to_bernstein(C, denominator, domain)
    biggest = max(abs(int(v)) for v in C.flat)
    precision = min(50, 58-(biggest.bit_length()-denominator.bit_length()+1))
    if precision < 0:
        return dict(status='UNKNOWN', reason='coefficient dynamic range')
    scale = 1 << precision
    low = np.array([int(v)*scale//denominator for v in C.flat], dtype=np.int64).reshape(C.shape)
    high = np.array([-((-int(v)*scale)//denominator) for v in C.flat], dtype=np.int64).reshape(C.shape)
    assert max(abs(low).max(), abs(high).max()) < 2**59
    del C
    angular_degree = core.gram()[2]

    def factor(point):
        value = 2*point[0]*point[1]
        if kind == 'full':
            value *= point[2]
        return value**angular_degree

    stack = [(low, high, ((0, 1),)*3, (0,)*3)]
    nodes = leaves = cap_leaves = maximum_depth = 0
    failures = []
    unknown = False
    while stack:
        lo, hi, box, levels = stack.pop()
        nodes += 1
        maximum_depth = max(maximum_depth, sum(levels))
        maximum = int(hi.max())
        if maximum <= 0:
            leaves += 1
            continue
        bounds = [(left+(right-left)*Q(x, 2**level), left+(right-left)*Q(y, 2**level))
                  for (left, right), (x, y), level in zip(domain, box, levels)]
        cap_lower, _ = correction_bounds(kind, bounds, regions, caps, core, basis_intervals)
        if Q(maximum, scale) <= factor([b[0] for b in bounds])*cap_lower:
            leaves += 1
            cap_leaves += 1
            continue
        failure = None
        for corner in itertools.product((0, 1), repeat=3):
            index = tuple(-1 if side else 0 for side in corner)
            base_lower = int(lo[index])
            if base_lower <= 0:
                continue
            point = [bounds[i][side] for i, side in enumerate(corner)]
            _, cap_upper = correction_bounds(kind, [(p, p) for p in point], regions,
                                             caps, core, basis_intervals)
            if Q(base_lower, scale) > factor(point)*cap_upper:
                failure = dict(point=list(map(str, point)), base_integer_lower=base_lower,
                               cap_upper=str(cap_upper))
                break
        if failure:
            failures.append(failure)
            break
        if nodes >= node_limit or sum(levels) >= depth_limit:
            unknown = True
            break
        axis = min(range(3), key=lambda j: levels[j])
        ll, lr = core.split(lo, axis, False)
        ul, ur = core.split(hi, axis, True)
        nextlevels = list(levels)
        nextlevels[axis] += 1
        leftbox, rightbox = list(box), list(box)
        a, b = box[axis]
        leftbox[axis], rightbox[axis] = (2*a, a+b), (a+b, 2*b)
        stack.extend([(lr, ur, tuple(rightbox), tuple(nextlevels)),
                      (ll, ul, tuple(leftbox), tuple(nextlevels))])
    status = 'VIOLATION' if failures else ('UNKNOWN' if unknown or stack else 'CERTIFIED')
    record = dict(status=status, kind=kind, regions=list(regions), nodes=nodes,
                  leaves=leaves, cap_leaves=cap_leaves, maxdepth=maximum_depth,
                  binary_scale=precision, failures=failures, seconds=time.time()-started,
                  pointed_helper='GENERAL_FULL_AND_EXTERIOR')
    print('EXACT FULL GENERAL POINTED', record, flush=True)
    return record
