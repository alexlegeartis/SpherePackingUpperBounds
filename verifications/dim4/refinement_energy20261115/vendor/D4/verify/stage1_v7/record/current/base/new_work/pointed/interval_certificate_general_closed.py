"""Extended exact pointed-ball enclosure for sigma in [-1,1].

The accepted v1 helper is unchanged. This module additionally supports probe
balls containing the packing origin and the R4 total-count bound five for
squared radius below 1/2. Its M is the active NONORIGIN count bound.
"""
from fractions import Fraction as Q
import itertools
import math
import time
import numpy as np
from interval_certificate import clipped, sqrt_bounds, radial_bounds


def probability_bounds(radius, cosine, cap):
    """Enclose the exact conditional ball probability on a rectangle."""
    lo, hi = radius
    tlo, thi = cosine
    d, a2, sigma = cap['d'], cap['a2'], cap['sigma']
    assert 0 < lo <= hi and d > 0 and a2 > 0
    assert -1 <= tlo <= thi <= 1 and -1 <= sigma <= 1
    b = d*d-a2
    c = lambda r: (r*r+b)/(2*d*r)
    if b <= 0:
        # c'(r)=(1-b/r^2)/(2d)>0.
        cmin, cmax = c(lo), c(hi)
    else:
        cmax = max(c(lo), c(hi))
        if lo*lo >= b:
            cmin = c(lo)
        elif hi*hi <= b:
            cmin = c(hi)
        else:
            cmin = sqrt_bounds(b)[0]/d
    assert cmin <= cmax
    if sigma == -1:
        if cmax+thi <= 0:
            return Q(1), Q(1)
        if cmin+tlo > 0:
            return Q(0), Q(0)
        return Q(0), Q(1)
    if sigma == 1:
        if tlo >= cmax:
            return Q(1), Q(1)
        if thi < cmin:
            return Q(0), Q(0)
        return Q(0), Q(1)
    products = [sigma*tlo, sigma*thi]
    Amin = cmin-max(products)
    Amax = cmax-min(products)
    minimum = min(1-tlo*tlo, 1-thi*thi)
    maximum = Q(1) if tlo <= 0 <= thi else max(1-tlo*tlo, 1-thi*thi)
    Dlo = sqrt_bounds((1-sigma*sigma)*minimum)[0]
    Dhi = sqrt_bounds((1-sigma*sigma)*maximum)[1]
    if Amax <= 0:
        lower = clipped((1-Amax/Dhi)/2) if Dhi else Q(1)
    else:
        lower = clipped((1-Amax/Dlo)/2) if Dlo else Q(0)
    # At A=D=0 the closed ball event has probability one, even when
    # the surrounding rectangle also has D>0.  The ordinary ratio formula
    # would otherwise give the false upper bound one half.
    if Amin == 0 and Dlo == 0:
        upper = Q(1)
    elif Amin >= 0:
        upper = clipped((1-Amin/Dhi)/2) if Dhi else (Q(0) if Amin > 0 else Q(1))
    else:
        upper = clipped((1-Amin/Dlo)/2) if Dlo else Q(1)
    assert 0 <= lower <= upper <= 1
    return lower, upper


def correction_bounds(bounds, regions, caps, core, basis_intervals):
    a, b, z = bounds
    cutoff = core.BREAKS[-1]
    assert 1 <= a[0] <= a[1] <= cutoff <= 2
    assert 1 <= b[0] <= b[1] <= cutoff and 0 <= z[0] <= z[1] <= 1
    h = lambda r, s: ((r+s)**2-cutoff**2)/(2*r*s)
    t = (-1+z[0]*h(a[0], b[0]), -1+z[1]*h(a[1], b[1]))
    assert -1 <= t[0] <= t[1] <= Q(1, 2)
    lower = upper = Q(0)
    for cap in caps:
        Pa = probability_bounds(a, t, cap)
        Pb = probability_bounds(b, t, cap)
        wa = radial_bounds(a, regions[0], cap, basis_intervals, core.DEG)
        wb = radial_bounds(b, regions[1], cap, basis_intervals, core.DEG)
        weight = cap['weight']/2
        lower += weight*(wb[0]*Pa[0]+wa[0]*Pb[0])
        upper += weight*(wb[1]*Pa[1]+wa[1]*Pb[1])
    return lower, upper


def validate_caps(records, core, basis_intervals):
    caps = []
    assert len(basis_intervals) == core.NS
    for i, (a, b) in enumerate(basis_intervals):
        assert a <= core.BREAKS[i] < core.BREAKS[i+1] <= b
    for record in records:
        cap = {key: Q(record[key]) for key in ['d', 'a2', 'sigma', 'weight']}
        cap['M'] = record['M']
        cap['w'] = list(map(Q, record['w']))
        assert isinstance(cap['M'], int) and cap['M'] >= 0
        assert cap['a2'] > 0 and cap['d'] > 0
        assert -1 <= cap['sigma'] <= 1 and cap['weight'] >= 0
        inside = int(cap['d']**2 <= cap['a2'])
        total = cap['M']+inside
        assert total >= 1
        if 'total_M' in record:
            assert isinstance(record['total_M'], int) and record['total_M'] == total
        if cap['a2'] < Q(total, 2*(total+1)):
            lemma = 'PAIR_DISTANCE_SUM'
        else:
            assert total >= 5 and cap['a2'] < Q(1, 2)
            lemma = 'R4_STRICT_NEGATIVE_INNER_PRODUCTS'
        cap.update(total_M=total, origin_inside=bool(inside), ball_count_lemma=lemma)
        assert len(cap['w']) == core.NC*core.NS and all(v >= 0 for v in cap['w'])
        caps.append(cap)
    return caps


def check(core, regions, eta, caps, basis_intervals, node_limit=180000, depth_limit=66):
    """Prove base Gamma minus generalized pointed corrections <=0."""
    start = time.time()
    C, denominator, domain = core.build('outside', regions, eta=eta)
    C, denominator = core.to_bernstein(C, denominator, domain)
    biggest = max(abs(int(v)) for v in C.flat)
    precision = min(50, 58-(biggest.bit_length()-denominator.bit_length()+1))
    if precision < 0:
        return {'status': 'UNKNOWN', 'reason': 'coefficient dynamic range'}
    scale = 1 << precision
    low = np.array([int(v)*scale//denominator for v in C.flat], dtype=np.int64).reshape(C.shape)
    high = np.array([-((-int(v)*scale)//denominator) for v in C.flat], dtype=np.int64).reshape(C.shape)
    assert max(abs(low).max(), abs(high).max()) < 2**59
    del C
    angular_degree = core.gram()[2]
    stack = [(low, high, ((0, 1),)*3, (0,)*3)]
    nodes = leaves = cap_leaves = maxdepth = 0
    failures = []
    unknown = False
    while stack:
        lo, hi, box, levels = stack.pop()
        nodes += 1
        maxdepth = max(maxdepth, sum(levels))
        maximum = int(hi.max())
        if maximum <= 0:
            leaves += 1
            continue
        bounds = [(left+(right-left)*Q(x, 2**level), left+(right-left)*Q(y, 2**level))
                  for (left, right), (x, y), level in zip(domain, box, levels)]
        cap_lower, _ = correction_bounds(bounds, regions, caps, core, basis_intervals)
        factor_lower = (2*bounds[0][0]*bounds[1][0])**angular_degree
        if Q(maximum, scale) <= factor_lower*cap_lower:
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
            _, cap_upper = correction_bounds([(p, p) for p in point], regions,
                                             caps, core, basis_intervals)
            if Q(base_lower, scale) > (2*point[0]*point[1])**angular_degree*cap_upper:
                failure = {'point': list(map(str, point)), 'base_integer_lower': base_lower,
                           'cap_upper': str(cap_upper)}
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
    report = {'status': status, 'regions': list(regions), 'nodes': nodes, 'leaves': leaves,
              'cap_leaves': cap_leaves, 'maxdepth': maxdepth, 'binary_scale': precision,
              'failures': failures, 'seconds': time.time()-start,
              'pointed_helper': 'GENERAL_SIGMA_AND_ACTIVE_BALL_COUNTS'}
    print('EXACT GENERAL POINTED EXTERIOR', report, flush=True)
    return report
