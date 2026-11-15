"""Rational Gram-ray and Bernstein certification for spherical S^3 bounds.

The triangle inequality is checked on the entire Gram-determinant domain.
Branches outside it are discarded only after an exact outward upper bound
for the determinant is strictly negative. UNKNOWN is never acceptance.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import argparse
import itertools
import json
import math
import hashlib
import time
import numpy as np


def convolve(a, b):
    result = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] += x*y
    return result


def radial_basis(degree, s, kind='bernstein'):
    p, q = s.numerator, s.denominator
    if kind == 'chebyshev':
        D = p+q
        numerator = [[1]]
        if degree:
            numerator.append([q-p, 2*q])
        for i in range(1, degree):
            next_row = [2*x for x in convolve([q-p, 2*q], numerator[-1])]
            for j, value in enumerate(numerator[-2]):
                next_row[j] -= D**2*value
            numerator.append(next_row)
        rows = [[x*D**(degree-i) for x in row]+[0]*(degree+1-len(row))
                for i, row in enumerate(numerator)]
        return np.array(rows, dtype=object), D**degree
    assert kind == 'bernstein'
    rows = []
    for i in range(degree+1):
        poly = [1]
        for _ in range(i):
            poly = convolve(poly, [q, q])
        for _ in range(degree-i):
            poly = convolve(poly, [p, -q])
        rows.append([math.comb(degree, i)*x for x in poly])
    return np.array(rows, dtype=object), (p+q)**degree


def gram_blocks(data):
    blocks = []
    nr = data['radial_degree']+1
    assert len(data['blocks']) == data['angular_degree']+1
    fixed_count = data.get('fixed_count', 0)
    if fixed_count:
        assert fixed_count in (23, 24, 25) and data['radial_basis'] == 'chebyshev'
        assert len(data['fixed_origin_free']) == nr+1
        assert all(isinstance(x, int) for x in data['fixed_origin_free'])
    for ell, rays in enumerate(data['blocks']):
        n = nr+1 if ell == 0 and not fixed_count else nr
        M = np.zeros((n, n), dtype=object)
        for vector in rays:
            assert len(vector) == n and all(isinstance(x, int) for x in vector)
            v = np.array(vector, dtype=object)
            M += np.outer(v, v)
        if ell == 0 and fixed_count:
            # M is the exact PSD matrix after substituting t_0=N-1. The
            # unrestricted physical origin block A can itself be indefinite.
            reduced = M
            m = fixed_count-1
            v0, b00, *h = data['fixed_origin_free']
            A = np.zeros((nr+1, nr+1), dtype=object)
            A[0, 0] = reduced[0, 0]-2*m*v0-m*m*b00
            A[0, 1] = A[1, 0] = v0
            A[1, 1] = b00
            A[1, 2:] = A[2:, 1] = h
            A[0, 2:] = A[2:, 0] = reduced[0, 1:]-m*np.array(h, dtype=object)
            A[2:, 2:] = reduced[1:, 1:]
            # Independent exact reconstruction of the restriction map.
            restriction = np.zeros((nr+1, nr), dtype=object)
            restriction[0, 0], restriction[1, 0] = 1, m
            for i in range(1, nr):
                restriction[i+1, i] = 1
            assert np.array_equal(restriction.T@A@restriction, reduced)
            M = A
        blocks.append(M)
    return blocks


def polynomials(data):
    assert data['dimension'] == 4
    s = F(data['s'])
    assert -1 < s < 1
    d, L = data['radial_degree'], data['angular_degree']
    assert len(data['scalar']) == data['scalar_degree']+1
    assert all(isinstance(x, int) for x in data['scalar']) and data['scalar_scale'] > 0
    assert data.get('verification_goal') == 'fixed_cardinality_signed_energy'
    assert data.get('fixed_count') in (23, 24, 25)
    # This goal is a bound only for the stated exact cardinality. The scalar
    # coefficients are unrestricted, and no all-cardinality claim is made.
    gram_den = data['feature_scale']**2
    assert gram_den > 0
    blocks = gram_blocks(data)
    radial, radial_den = radial_basis(d, s, data.get('radial_basis', 'bernstein'))
    powers = [radial.T@(M[1:, 1:] if ell == 0 else M)@radial
              for ell, M in enumerate(blocks)]
    size = d+L+1
    triangle = np.zeros((size, size, size), dtype=object)
    triangle_den = gram_den*radial_den**2*2**L
    for ell, M in enumerate(powers):
        for j in range(ell//2+1):
            angular = ((-1)**j*math.comb(ell, j)*math.comb(2*ell-2*j, ell)*2**(L-ell))
            k = ell-2*j
            for h in range(k+1):
                for u in range(j+1):
                    for v in range(j+1):
                        coef = angular*math.comb(k, h)*math.comb(j, u)*math.comb(j, v)*(-1)**(h+u+v)
                        ea, eb, et = h+2*u, h+2*v, k-h
                        for i in range(d+1):
                            for q in range(d+1):
                                term = coef*M[i, q]
                                ia, ib = i+ea, q+eb
                                triangle[ia, ib, et] += term
                                triangle[ia, et, ib] += term
                                triangle[et, ia, ib] += term
    assert np.array_equal(triangle, np.swapaxes(triangle, 0, 1))
    assert np.array_equal(triangle, np.swapaxes(triangle, 1, 2))
    # Tapered harmonic blocks can make the actual degree far smaller than
    # d+L. Remove only identically zero coefficient faces, exactly.
    while triangle.shape[0] > 1 and not np.any(triangle[-1, :, :]):
        triangle = triangle[:-1, :-1, :-1]
    maxdegree = max(data['scalar_degree'], 2*d+2*L)
    pair = [F(0)]*(maxdegree+1)
    for ell, a in enumerate(data['scalar']):
        for j in range(ell//2+1):
            power = ell-2*j
            pair[power] += F(a*(-1)**j*math.comb(ell-j, j)*2**power,
                             data['scalar_scale']*(ell+1))
    coupling = blocks[0][0, 1:]@radial
    for i, x in enumerate(coupling):
        pair[i] += F(2*int(x), gram_den*radial_den)
    for ell, M in enumerate(powers):
        for i in range(d+1):
            for j in range(d+1):
                for k in range(ell+1):
                    pair[i+j+2*k] += F(int(M[i, j])*math.comb(ell, k)*(-1)**k,
                                       gram_den*radial_den**2)
    pair_den = math.lcm(*(x.denominator for x in pair))
    while len(pair) > 1 and pair[-1] == 0:
        pair.pop()
    pair = np.array([int(x*pair_den) for x in pair], dtype=object)
    objective = F(sum(data['scalar']), data['scalar_scale'])+F(int(blocks[0][0, 0]), gram_den)
    determinant = np.zeros((3, 3, 3), dtype=object)
    determinant[0, 0, 0] = 1
    determinant[1, 1, 1] = 2
    determinant[2, 0, 0] = determinant[0, 2, 0] = determinant[0, 0, 2] = -1
    return triangle, triangle_den, pair, pair_den, determinant, objective, s


def audit_identities(data, triangle, td, pair, pd):
    """Exact independent direct-kernel evaluations at rational arguments."""
    d, L, s = data['radial_degree'], data['angular_degree'], F(data['s'])
    blocks = gram_blocks(data)
    den = data['feature_scale']**2

    def b(a):
        z = (a+1)/(s+1)
        if data.get('radial_basis', 'bernstein') == 'chebyshev':
            out = [F(1)]
            if d:
                out.append(2*z-1)
            for i in range(1, d):
                out.append(2*(2*z-1)*out[-1]-out[-2])
            return out
        return [math.comb(d, i)*z**i*(1-z)**(d-i) for i in range(d+1)]

    def kernel(a, b_, t):
        ba, bb = b(a), b(b_)
        x, w = t-a*b_, (1-a*a)*(1-b_*b_)
        H = [F(1)]
        if L:
            H.append(x)
        for ell in range(1, L):
            H.append(((2*ell+1)*x*H[-1]-ell*w*H[-2])/(ell+1))
        answer = F(0)
        for ell, M in enumerate(blocks):
            M = M[1:, 1:] if ell == 0 else M
            answer += H[ell]*sum((F(int(M[i, j]), den)*ba[i]*bb[j]
                                  for i in range(d+1) for j in range(d+1)), F(0))
        return answer

    def tensor_value(C, point):
        result = C
        denominator = 1
        weights = []
        for axis, x in enumerate(point):
            n = C.shape[axis]-1
            weights.append(np.array([x.numerator**j*x.denominator**(n-j)
                                     for j in range(n+1)], dtype=object))
            denominator *= x.denominator**n
        for weight in reversed(weights):
            result = np.tensordot(result, weight, axes=(-1, 0))
        return F(int(result), denominator)

    for a, b_, t in [(F(-1), F(-1, 2), F(1, 2)),
                      (F(-1, 3), F(1, 5), F(2, 5))]:
        direct = kernel(a, b_, t)+kernel(a, t, b_)+kernel(b_, t, a)
        assert tensor_value(triangle, [a, b_, t])/td == direct
    for a in [F(-1), F(0), s]:
        U = [F(1)]
        if data['scalar_degree']:
            U.append(2*a)
        for ell in range(2, data['scalar_degree']+1):
            U.append(2*a*U[-1]-U[-2])
        scalar = sum(F(x, data['scalar_scale'])*U[j]/(j+1)
                     for j, x in enumerate(data['scalar']))
        ba = b(a)
        coupling = sum(F(int(blocks[0][0, i+1]), den)*ba[i] for i in range(d+1))
        direct = scalar+2*coupling+kernel(a, a, F(1))
        assert tensor_value(pair, [a])/pd == direct
    print('PASS five exact direct-kernel/power-polynomial identities.', flush=True)


@lru_cache(None)
def transform(degree, left, width):
    rows = []
    for i in range(degree+1):
        row = []
        for k in range(degree+1):
            row.append(sum((F(math.comb(k, j)*math.comb(i, j), math.comb(degree, j))
                            *left**(k-j)*width**j for j in range(min(i, k)+1)), F(0)))
        rows.append(row)
    denominator = math.lcm(*(x.denominator for row in rows for x in row))
    return np.array([[int(x*denominator) for x in row] for row in rows], dtype=object), denominator


def bernstein(coefficients, denominator, bounds):
    coefficients = coefficients.copy()
    for axis, (left, right) in enumerate(bounds):
        n = coefficients.shape[axis]
        T, td = transform(n-1, left, right-left)
        moved = np.moveaxis(coefficients, axis, 0)
        shape = moved.shape
        coefficients = np.moveaxis((T@moved.reshape(n, -1)).reshape(shape), 0, axis)
        denominator *= td
    big = max(abs(int(x)) for x in coefficients.flat)
    precision = min(48, 58-(big.bit_length()-denominator.bit_length()+1))
    assert precision >= 0, 'Coefficient range exceeds fixed-point budget'
    lower = np.array([int(x)*2**precision//denominator for x in coefficients.flat],
                     dtype=np.int64).reshape(coefficients.shape)
    upper = np.array([-((-int(x)*2**precision)//denominator) for x in coefficients.flat],
                     dtype=np.int64).reshape(coefficients.shape)
    assert max(abs(lower).max(), abs(upper).max()) < 2**59
    return lower, upper, precision


def split(coefficients, axis, upper):
    work = np.moveaxis(coefficients, axis, 0).copy()
    n = work.shape[0]-1
    left, right = np.empty_like(work), np.empty_like(work)
    left[0], right[n] = work[0], work[n]
    for j in range(1, n+1):
        work = (work[:-1]+work[1:]+int(upper))//2
        left[j], right[n-j] = work[0], work[-1]
    return np.moveaxis(left, 0, axis), np.moveaxis(right, 0, axis)


def certify(coefficients, denominator, bounds, gram=None, node_limit=250000, depth_limit=66):
    low, high, precision = bernstein(coefficients, denominator, bounds)
    if gram is None:
        glo = ghi = gp = None
    else:
        glo, ghi, gp = bernstein(gram, 1, bounds)
    dims = len(bounds)
    stack = [(low, high, glo, ghi, (0,)*dims, ((0, 1),)*dims)]
    nodes = accepted = excluded = depth = 0
    while stack:
        lo, hi, gl, gh, levels, box = stack.pop()
        nodes += 1
        depth = max(depth, sum(levels))
        if gh is not None and gh.max() < 0:
            excluded += 1
            continue
        if hi.max() <= 0:
            accepted += 1
            continue
        for corner in itertools.product((0, 1), repeat=dims):
            index = tuple(-1 if x else 0 for x in corner)
            if lo[index] > 0 and (gl is None or gl[index] >= 0):
                point = [bounds[i][0]+(bounds[i][1]-bounds[i][0])*
                         F(box[i][corner[i]], 2**levels[i]) for i in range(dims)]
                return {'status': 'VIOLATION', 'point': list(map(str, point)),
                        'nodes': nodes, 'depth': depth, 'lower_integer': int(lo[index]),
                        'binary_precision': precision}
        if nodes >= node_limit or sum(levels) >= depth_limit:
            return {'status': 'UNKNOWN', 'nodes': nodes, 'depth': depth,
                    'accepted': accepted, 'excluded': excluded,
                    'box': box, 'levels': levels}
        axis = min(range(dims), key=lambda i: levels[i])
        ll, lr = split(lo, axis, False)
        ul, ur = split(hi, axis, True)
        if gl is not None:
            gll, glr = split(gl, axis, False)
            gul, gur = split(gh, axis, True)
        else:
            gll = glr = gul = gur = None
        next_levels = list(levels)
        next_levels[axis] += 1
        next_levels = tuple(next_levels)
        bl, br = list(box), list(box)
        a, b = box[axis]
        bl[axis], br[axis] = (2*a, a+b), (a+b, 2*b)
        stack.extend([(lr, ur, glr, gur, next_levels, tuple(br)),
                      (ll, ul, gll, gul, next_levels, tuple(bl))])
    return {'status': 'CERTIFIED', 'nodes': nodes, 'accepted': accepted,
            'excluded': excluded, 'maxdepth': depth, 'binary_precision': precision}
