"""Exact pair majorant checks retaining negative pointed self diagonals.

The caller supplies an exact polynomial C/den majorizing
    F(alpha*r^2) + exp(alpha*r^2)*D_restored(r),
where D_restored includes the PSD diagonal, restored origin row, and count
terms.  This checker proves C/den <= exp(alpha*r^2)*J(r), where J is the
nonnegative magnitude of the pointed self diagonal.  Together these two
inequalities imply the actual pair condition on the entire interval.

Failure means this sufficient comparison was not proved; a positive
majorant witness alone is not labelled a violation of the original pair
condition.  Accepted v1 source files are left unchanged.
"""
from fractions import Fraction as Q
import math
import time
import numpy as np
from interval_certificate import radial_bounds
from interval_certificate_general_closed import probability_bounds


def exponential_bounds(t, order=55):
    assert 0 <= t < order+2
    partial = Q(1)
    term = Q(1)
    for k in range(1, order+1):
        term *= t/k
        partial += term
    nextterm = term*t/(order+1)
    # Successive omitted terms have ratio at most t/(order+2).
    return partial, partial+nextterm/(1-t/(order+2))


def diagonal_bounds(radius, region, caps, core, basis_intervals):
    """Enclose the nonnegative reduction from the auxiliary pair row.

    An untruncated ray reduces by its self diagonal.  A ray truncated below
    C reduces by M*w at r<=C, and by its usual self diagonal at r>C.  This
    represents the actual truncated linear coupling exactly.
    """
    lower = upper = Q(0)
    for cap in caps:
        p = probability_bounds(radius, (Q(1), Q(1)), cap)
        if 'pivot_lower_cutoff' in cap:
            cutoff = cap['pivot_lower_cutoff']
            if radius[1] <= cutoff:
                p = (Q(cap['M']), Q(cap['M']))
            elif radius[0] <= cutoff:
                p = (min(Q(cap['M']), p[0]), max(Q(cap['M']), p[1]))
        if p[1] == 0:
            continue
        w = radial_bounds(radius, region, cap, basis_intervals, core.DEG)
        lower += cap['weight']*p[0]*w[0]
        upper += cap['weight']*p[1]*w[1]
    return lower, upper


def check(core, polynomial, denominator, interval, region, caps,
          basis_intervals, alpha=Q(2), node_limit=120000, depth_limit=66):
    started = time.time()
    assert alpha > 0 and 0 < interval[0] <= interval[1]
    C, denominator = core.to_bernstein(polynomial, denominator, [interval])
    biggest = max(abs(int(v)) for v in C.flat)
    precision = min(50, 58-(biggest.bit_length()-denominator.bit_length()+1))
    if precision < 0:
        return dict(status='UNKNOWN', reason='pair coefficient dynamic range')
    scale = 1 << precision
    low = np.array([int(v)*scale//denominator for v in C.flat], dtype=np.int64)
    high = np.array([-((-int(v)*scale)//denominator) for v in C.flat], dtype=np.int64)
    assert max(abs(low).max(), abs(high).max()) < 2**59
    del C
    stack = [(low, high, interval, 0)]
    nodes = leaves = diagonal_leaves = maximum_depth = 0
    failures = []
    unknown = False
    while stack:
        lo, hi, radius, depth = stack.pop()
        nodes += 1
        maximum_depth = max(maximum_depth, depth)
        maximum = int(hi.max())
        if maximum <= 0:
            leaves += 1
            continue
        diagonal_lower, _ = diagonal_bounds(radius, region, caps, core, basis_intervals)
        if diagonal_lower:
            exp_lower = exponential_bounds(alpha*radius[0]**2)[0]
            if Q(maximum, scale) <= exp_lower*diagonal_lower:
                leaves += 1
                diagonal_leaves += 1
                continue
        failure = None
        for i, r in [(0, radius[0]), (-1, radius[1])]:
            if lo[i] <= 0:
                continue
            _, diagonal_upper = diagonal_bounds((r, r), region, caps, core, basis_intervals)
            exp_upper = exponential_bounds(alpha*r*r)[1]
            if Q(int(lo[i]), scale) > exp_upper*diagonal_upper:
                failure = dict(radius=str(r), integer_majorant_lower=int(lo[i]),
                               exponential_upper=str(exp_upper),
                               pointed_diagonal_upper=str(diagonal_upper),
                               failure_kind='SUFFICIENT_PAIR_MAJORANT_COMPARISON')
                break
        if failure:
            failures.append(failure)
            break
        if nodes >= node_limit or depth >= depth_limit:
            unknown = True
            break
        ll, lr = core.split(lo, 0, False)
        ul, ur = core.split(hi, 0, True)
        midpoint = sum(radius)/2
        stack.extend([(lr, ur, (midpoint, radius[1]), depth+1),
                      (ll, ul, (radius[0], midpoint), depth+1)])
    status = 'VIOLATION' if failures else ('UNKNOWN' if unknown or stack else 'CERTIFIED')
    result = dict(status=status, region=region, left=str(interval[0]),
                  right=str(interval[1]), nodes=nodes, leaves=leaves,
                  diagonal_leaves=diagonal_leaves, maxdepth=maximum_depth,
                  binary_scale=precision, failures=failures,
                  seconds=time.time()-started,
                  pointed_pair_contract='EXACT_SELF_DIAGONAL_AND_VALIDATED_PIVOT_TRUNCATION')
    print('EXACT POINTED PAIR MAJORANT', result, flush=True)
    return result
