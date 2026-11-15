"""Exact cap compression by affine latitude, checked on the full domain.

The spherical transformation fixes the south pole and sends latitude a to
phi(a)=alpha*a-beta, alpha=1-beta.  For fixed original inner product s,
its transformed inner product is phi(a)phi(b)+sqrt(K(a)K(b))*(s-a*b),
where K(a)=alpha*((1-a)+beta*(1+a))/(1-a).

The radical comparison is reduced to rational polynomials P,V,H and checked
on every feasible latitude pair using integer Bernstein subdivision.
The numerical optimizer and any conjectured extremizer play no role here.
This file establishes conditional local counts, not a global density bound.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse
import hashlib
import itertools
import json
import math
import time

HERE = Path(__file__).resolve().parent


def add(*polynomials):
    out = {}
    for p in polynomials:
        for key, value in p.items():
            out[key] = out.get(key, Q(0))+value
    return {key:value for key, value in out.items() if value}


def mul(a, b):
    out = {}
    for (i, j), x in a.items():
        for (k, ell), y in b.items():
            key = (i+k, j+ell)
            out[key] = out.get(key, Q(0))+x*y
    return {key:value for key, value in out.items() if value}


def scale(a, value):
    return {key:x*value for key, x in a.items() if x*value}


def polynomials(s, c, S):
    beta = max(Q(0), (c-S)/(1+c)); alpha = 1-beta
    assert 0 <= beta < 1 and -1 < s < 1 and c < 1
    a, b, one = {(1, 0):Q(1)}, {(0, 1):Q(1)}, {(0, 0):Q(1)}
    ab = mul(a, b)
    phi_a, phi_b = add(scale(a, alpha), scale(one, -beta)), add(scale(b, alpha), scale(one, -beta))
    P = add(scale(one, S), scale(mul(phi_a, phi_b), -1))
    V = add(scale(one, s), scale(ab, -1))
    ma, mb = add(one, scale(a, -1)), add(one, scale(b, -1))
    ka = add(ma, scale(add(one, a), beta))
    kb = add(mb, scale(add(one, b), beta))
    H = add(mul(mul(P, P), mul(ma, mb)),
            scale(mul(mul(ka, kb), mul(V, V)), -alpha*alpha))
    D = add(scale(one, 1-s*s), scale(ab, 2*s), scale(mul(a, a), -1), scale(mul(b, b), -1))
    return dict(P=P, V=V, H=H, D=D), beta


def bernstein(poly, lo, hi):
    degree = max(max(key) for key in poly)
    width = hi-lo
    power = [[Q(0)]*(degree+1) for _ in range(degree+1)]
    for (i, j), value in poly.items():
        for k in range(i+1):
            for ell in range(j+1):
                power[k][ell] += value*math.comb(i, k)*math.comb(j, ell)*lo**(i+j-k-ell)*width**(k+ell)
    values = [[sum((power[k][ell]*Q(math.comb(i, k), math.comb(degree, k))*
                             Q(math.comb(j, ell), math.comb(degree, ell))
                             for k in range(i+1) for ell in range(j+1)), Q(0))
               for j in range(degree+1)] for i in range(degree+1)]
    denominator = math.lcm(*(v.denominator for row in values for v in row))
    return [[int(v*denominator) for v in row] for row in values]


def split(coeff, axis):
    n = len(coeff); d = n-1
    left, right = [[0]*n for _ in range(n)], [[0]*n for _ in range(n)]
    for line in range(n):
        values = [coeff[i][line] if axis == 0 else coeff[line][i] for i in range(n)]
        for i in range(n):
            lo = (1 << (d-i))*sum(math.comb(i, j)*values[j] for j in range(i+1))
            hi = (1 << i)*sum(math.comb(d-i, j-i)*values[j] for j in range(i, n))
            if axis == 0:
                left[i][line], right[i][line] = lo, hi
            else:
                left[line][i], right[line][i] = lo, hi
    return left, right


def bounds(coeff):
    return min(map(min, coeff)), max(map(max, coeff))


def check(inner, outer, S=Q(203, 400), max_nodes=2000000, max_depth=100):
    if not __debug__:
        raise RuntimeError('Run without -O.')
    inner, outer, S = Q(inner), Q(outer), Q(S)
    assert 1 <= inner < outer and outer*outer < 2
    s = 1-1/(2*inner*inner)
    c = (inner*inner+outer*outer-1)/(2*inner*outer)
    polys, beta = polynomials(s, c, S)
    initial = {key:bernstein(p, Q(-1), c) for key, p in polys.items()}
    queue = [(initial, 0)]; nodes = 0; accepted = [0, 0, 0]; excluded = 0; depth_used = 0
    started = time.time()
    while queue:
        coeff, depth = queue.pop(); nodes += 1; depth_used = max(depth_used, depth)
        if nodes > max_nodes or depth > max_depth:
            return dict(status='UNKNOWN', nodes=nodes, maxdepth=depth_used,
                        seconds=time.time()-started, inner_radius=str(inner), outer_radius=str(outer))
        loD, hiD = bounds(coeff['D'])
        if hiD < 0:
            excluded += 1
            continue
        loP, hiP = bounds(coeff['P']); loV, hiV = bounds(coeff['V']); loH, hiH = bounds(coeff['H'])
        if loP >= 0 and hiV <= 0:
            accepted[0] += 1
            continue
        if loP >= 0 and loH >= 0:
            accepted[1] += 1
            continue
        if hiV <= 0 and hiH <= 0:
            accepted[2] += 1
            continue
        children = [{}, {}]
        for key, value in coeff.items():
            halves = split(value, depth % 2)
            children[0][key], children[1][key] = halves
        queue.extend((child, depth+1) for child in children)
    return dict(status='EXACT_AFFINE_LATITUDE_DOMAIN_CERTIFIED',
                inner_radius=str(inner), outer_radius=str(outer),
                spherical_cutoff=str(S), inner_cosine_bound=str(s), pole_cosine_bound=str(c),
                beta=str(beta), alpha=str(1-beta),
                transformed_pole_upper=str((1-beta)*c-beta),
                nodes=nodes, accepted_sign_cases=accepted, excluded_infeasible=excluded,
                maxdepth=depth_used, seconds=time.time()-started,
                polynomial_degrees={key:max(max(v) for v in p) for key, p in polys.items()},
                sign_cases=['P>=0,V<=0', 'P>=0,H>=0', 'V<=0,H<=0'])


def verify_record(record, verified_spherical):
    assert verified_spherical['status'] == 'EXACT_SPHERICAL_COUNT_24_VERIFIED'
    assert verified_spherical['count_bound'] <= 24
    assert Q(record['spherical_cutoff']) <= Q(verified_spherical['cosine_cutoff'])
    checked = check(record['inner_radius'], record['outer_radius'], Q(record['spherical_cutoff']))
    assert checked['status'] == 'EXACT_AFFINE_LATITUDE_DOMAIN_CERTIFIED'
    for key in ['inner_cosine_bound', 'pole_cosine_bound', 'beta', 'alpha', 'transformed_pole_upper']:
        assert checked[key] == record[key]
    assert Q(checked['transformed_pole_upper']) <= Q(record['spherical_cutoff'])
    return checked


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inner', type=Q)
    parser.add_argument('--outer', type=Q)
    parser.add_argument('--max-nodes', type=int, default=2000000)
    parser.add_argument('--skip-spherical-replay', action='store_true')
    args = parser.parse_args()
    pairs = [(Q('1.005643'), Q('1.01534')), (Q('1.002946'), Q('1.02597'))]
    if args.inner is not None:
        assert args.outer is not None
        pairs = [(args.inner, args.outer)]
    reports = []
    for inner, outer in pairs:
        report = check(inner, outer, max_nodes=args.max_nodes)
        reports.append(report)
        print(json.dumps(report), flush=True)
    result = dict(status='EXACT_AFFINE_CONDITIONAL_GEOMETRY_VERIFIED' if all(
                      r['status']=='EXACT_AFFINE_LATITUDE_DOMAIN_CERTIFIED' for r in reports) else 'INCOMPLETE',
                  spherical_dependency='A(4,203/400)<=24; separately retained exact certificate',
                  records=reports, verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  note='Only the full latitude comparison is verified here; the spherical proof is a separate dependency.')
    if not args.skip_spherical_replay and result['status'] != 'INCOMPLETE':
        from mobius_counts import replay_spherical, CERTIFICATE
        spherical = replay_spherical()
        result['spherical_verification'] = spherical
        result['spherical_certificate'] = CERTIFICATE.name
        result['spherical_certificate_sha256'] = hashlib.sha256(CERTIFICATE.read_bytes()).hexdigest()
        assert spherical['certificate_sha256'] == result['spherical_certificate_sha256']
        for record in reports:
            assert Q(record['spherical_cutoff']) <= Q(spherical['cosine_cutoff'])
            assert spherical['status'] == 'EXACT_SPHERICAL_COUNT_24_VERIFIED' and spherical['count_bound'] <= 24
            assert Q(record['transformed_pole_upper']) <= Q(record['spherical_cutoff'])
            record['maximum_inner_count'] = 24
            record['spherical_certificate'] = CERTIFICATE.name
            record['spherical_certificate_sha256'] = result['spherical_certificate_sha256']
            record['local_energy'] = '(23-N_inner)*(N_outer-N_inner) >= 0'
        result['status'] = 'EXACT_AFFINE_CONDITIONAL_COUNTS_VERIFIED_WITH_SPHERICAL_REPLAY'
        result['note'] = 'Exact local count lemmas only; a separate global certificate is required for a density bound.'
    target = HERE/'affine_latitude_certificates.json'
    target.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], flush=True)


if __name__ == '__main__':
    main()
