"""Exact cyclic verification with a previously proved full pair/count gauge.

This separate checker extends the audited closed-boundary cap/tent/count
helper by subtracting the independently validated antipodal correction.
The allocation is also present in the complete pair budget, via the
mandatory linear-origin identity in its calling global verifier.
The positive Gaussian-polynomial budget must first be established by the
calling verifier's fresh exact unlifted pair checks.  Accepted historical
helper sources are not modified.
"""
from fractions import Fraction as Q
import itertools
import time
import numpy as np
from interval_certificate_tent_closed import correction_bounds as tent_bounds
from interval_certificate_full_closed import correction_bounds as cap_bounds
from full_pair_count_gauge import correction_bounds as lift_bounds
from two_capacity_antipodal_allocation_exact import correction_bounds as allocation_bounds


def check(core, kind, regions, eta, tents, basis_intervals, budget, probes,
          node_limit=180000, depth_limit=66, caps=(), allocations=()):
    """Prove one closed cyclic branch, with no inference about other branches."""
    started = time.time()
    assert kind in ['full', 'outside'] and budget.core is core
    assert budget.nonnegative_proven, 'Fresh exact base-pair proofs are required first.'
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

    def all_corrections(bounds):
        tl, tu = tent_bounds(kind, bounds, regions, tents, core, basis_intervals)
        cl, cu = cap_bounds(kind, bounds, regions, caps, core, basis_intervals) if caps else (Q(0), Q(0))
        ll, lu = lift_bounds(kind, bounds, regions, probes, budget)
        al, au = allocation_bounds(kind, bounds, regions, allocations, core,
                                   basis_intervals)
        assert min(tl, cl, ll, al) >= 0
        return tl+cl+ll+al, tu+cu+lu+au, tl+cl+al, tl+cl+ll

    stack = [(low, high, ((0, 1),)*3, (0,)*3)]
    nodes = leaves = correction_leaves = lift_dependent_leaves = maximum_depth = 0
    antipodal_dependent_leaves = 0
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
        correction_lower, _, unlifted_lower, without_antipodal = all_corrections(bounds)
        factor_lower = factor([b[0] for b in bounds])
        if Q(maximum, scale) <= factor_lower*correction_lower:
            leaves += 1
            correction_leaves += 1
            lift_dependent_leaves += int(Q(maximum, scale) > factor_lower*unlifted_lower)
            antipodal_dependent_leaves += int(Q(maximum, scale) > factor_lower*without_antipodal)
            continue
        failure = None
        for corner in itertools.product((0, 1), repeat=3):
            index = tuple(-1 if side else 0 for side in corner)
            base_lower = int(lo[index])
            if base_lower <= 0:
                continue
            point = [bounds[i][side] for i, side in enumerate(corner)]
            _, correction_upper, _, _ = all_corrections([(p, p) for p in point])
            if Q(base_lower, scale) > factor(point)*correction_upper:
                failure = dict(point=list(map(str, point)), base_integer_lower=base_lower,
                               correction_upper=str(correction_upper))
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
    record = dict(status=status, kind=kind, regions=list(regions), nodes=nodes, leaves=leaves,
                  correction_leaves=correction_leaves, lift_dependent_leaves=lift_dependent_leaves,
                  antipodal_dependent_leaves=antipodal_dependent_leaves,
                  maxdepth=maximum_depth, binary_scale=precision, failures=failures,
                  seconds=time.time()-started,
                  pointed_helper='EXACT_CAP_TENT_TWO_CAPACITY_ANTIPODAL_AND_VERIFIED_FULL_PAIR_COUNT_GAUGE')
    print('EXACT TWO CAPACITY ANTIPODAL COUNT GAUGE CYCLIC', record, flush=True)
    return record
