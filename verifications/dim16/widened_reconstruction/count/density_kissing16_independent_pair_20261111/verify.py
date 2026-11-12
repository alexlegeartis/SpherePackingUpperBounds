#!/usr/bin/env python3
"""Independent exact rational PSD, origin and pair verifier.

This script imports neither discovery code nor any previous verifier. Every JSON
coefficient is parsed from its printed decimal representation, without floats.
The cyclic/three-point gate is deliberately outside its scope.
"""
if not __debug__:
    raise RuntimeError("Certificate verification requires assertions; do not run Python with -O or -OO.")

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from time import monotonic
from flint import fmpq, fmpq_poly, ctx
ctx.threads = 1

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
FROZEN = HERE.parent / 'density_kissing16_reweight_frozen_20261111'
EXPECTED = {
    'candidate.json': 'ec45e2946fc318362649b92d519c915f26c73e524fc284b1ece0940e8bcdd257',
    'factors.json': 'a190925c6f7e781b5905064f9aebfdfc938d8a25a82a27ace93e7ea74aafd062',
}

def sha(data):
    return hashlib.sha256(data).hexdigest()

def Q(value):
    f = Fraction(str(value))
    return fmpq(f.numerator, f.denominator)

def dump(name, value):
    (HERE / name).write_text(json.dumps(value, indent=2) + '\n')

def decimal_enclosure(x, digits=24):
    n, d = int(x.numer()), int(x.denom())
    scale = 10**digits
    lower = (n * scale) // d
    upper = -((-n * scale) // d)
    def fmt(n):
        sign = '-' if n < 0 else ''
        s = str(abs(n)).rjust(digits + 1, '0')
        return sign + s[:-digits] + '.' + s[-digits:]
    return {'lower': fmt(lower), 'upper': fmt(upper), 'digits': digits}

def polys(n, degree):
    # P_0=1, P_1=x; normalized Gegenbauer recurrence with P_j(1)=1.
    p = [fmpq_poly([1])]
    if degree:
        p.append(fmpq_poly([0, 1]))
    x = fmpq_poly([0, 1])
    for j in range(1, degree):
        p.append(((2*j+n-2)*x*p[-1] - j*p[-2]) / (j+n-2))
    assert all(a(1) == 1 for a in p)
    assert all(a(-1) == (-1)**j for j, a in enumerate(p))
    return p

def psd_blocks(factors):
    blocks = factors['blocks']
    assert len(blocks) == 34
    parsed = []
    records = 0
    for k, b in enumerate(blocks):
        assert b['index'] == k
        size = 1 if k < 16 else (2 if k == 16 else 34-k)
        assert b['size'] == size
        lam = Q(b['lambda_decimal'])
        assert lam > 0
        L = [[Q(0) for _ in range(size)] for _ in range(size)]
        seen = set()
        for i, j, value in b['L_entries']:
            assert 0 <= i < size and 0 <= j < size and (i, j) not in seen
            assert isinstance(value, str)
            seen.add((i,j))
            L[i][j] = Q(value)
            records += 1
        # Independent dense reconstruction of X = L^T L + lambda I.
        X = [[sum((L[r][i]*L[r][j] for r in range(size)), Q(0))
              + (lam if i == j else Q(0))
              for j in range(size)] for i in range(size)]
        assert all(X[i][j] == X[j][i] for i in range(size) for j in range(size))
        parsed.append((L, lam, X))
    assert records == factors['factor_records']
    return parsed, records

def factor_quadratic(L, lam, p):
    transformed = [sum((v*q for v, q in zip(row, p)), fmpq_poly([])) for row in L]
    return sum((q*q for q in transformed), fmpq_poly([])) + lam*sum((q*q for q in p), fmpq_poly([]))

def dense_quadratic(X, p):
    return sum((X[i][j]*p[i]*p[j] for i in range(len(p)) for j in range(len(p))), fmpq_poly([]))

def factor_cross(L, lam, p):
    return sum((sum(row)*sum((v*q for v,q in zip(row,p)), fmpq_poly([])) for row in L), fmpq_poly([])) + lam*sum(p, fmpq_poly([]))

def primitive_positive_scale(p):
    """Remove only positive scalar factors, preserving every pointwise sign."""
    n = p.numer()
    if n.is_zero():
        return p
    content = abs(n.content())
    assert content > 0 and p.denom() > 0
    return fmpq_poly([fmpq(v, content) for v in n.coeffs()])

def sign(x):
    return int(x > 0) - int(x < 0)

def variations(signs):
    nz = [s for s in signs if s]
    return sum(a != b for a,b in zip(nz,nz[1:]))

def sturm_certificate(p, left, right, name):
    """Exact sign-preserving Sturm sequence; root-free + negative proves < 0."""
    chain = [primitive_positive_scale(p), primitive_positive_scale(p.derivative())]
    while not chain[-1].is_constant():
        remainder = -(chain[-2] % chain[-1])
        if remainder.is_zero():
            break
        chain.append(primitive_positive_scale(remainder))
    signs_left = [sign(q(left)) for q in chain]
    signs_right = [sign(q(right)) for q in chain]
    va, vb = variations(signs_left), variations(signs_right)
    assert p(left) < 0 and p(right) < 0
    assert va-vb == 0
    assert chain[-1].degree() == 0  # Squarefree; no unreported endpoint convention.
    # Retain all integer primitive polynomials, allowing simple independent replay.
    serial = [[str(c) for c in q.coeffs()] for q in chain]
    raw = (json.dumps(serial, separators=(',', ':')) + '\n').encode()
    certificate_file = name + '_sturm_chain.json'
    (HERE / certificate_file).write_bytes(raw)
    return {
        'polynomial_shift': str(p - PAIR),
        'degree': p.degree(),
        'chain_degrees': [q.degree() for q in chain],
        'left_signs': signs_left,
        'right_signs': signs_right,
        'left_variations': va,
        'right_variations': vb,
        'distinct_roots_in_open_interval': va-vb,
        'left_exact': str(p(left)),
        'right_exact': str(p(right)),
        'left_decimal_enclosure': decimal_enclosure(p(left)),
        'right_decimal_enclosure': decimal_enclosure(p(right)),
        'chain_file': certificate_file,
        'chain_sha256': sha(raw),
        'conclusion': 'Strictly negative on the entire closed interval by exact Sturm theorem and continuity.',
    }

def main():
    global PAIR
    start = monotonic()
    inputs = {}
    for name, expected in EXPECTED.items():
        data = (FROZEN / name).read_bytes()
        assert sha(data) == expected, (name, sha(data))
        inputs[name] = json.loads(data, parse_float=str)
    candidate, factors = inputs['candidate.json'], inputs['factors.json']
    assert candidate['factor_sha256'] == EXPECTED['factors.json']
    assert candidate['dimension'] == factors['dimension'] == 16
    assert factors['degree_parameter'] == 16
    assert candidate['cap_exact'] == '501/1000'
    left, cap = Q(-1), Q(candidate['cap_exact'])
    blocks, records = psd_blocks(factors)
    rec = candidate['record']
    theta = list(map(Q, rec['theta']))
    a = list(map(Q, rec['global_coefficients']))
    r = list(map(Q, rec['local_angular_coefficients']))
    assert len(theta) == 17 and len(a) == len(r) == 24
    assert all(v >= 0 for v in theta+a+r)
    assert theta[0] > 0
    B = blocks[16][2]
    p_global = polys(16,24)
    f = fmpq_poly([1]) + sum((a[j-1]*p_global[j] for j in range(1,25)), fmpq_poly([]))
    R = sum((r[j-1]*p_global[j] for j in range(1,25)), fmpq_poly([]))
    assert f(1) == 1+sum(a) and R(1) == sum(r)
    # Homogeneous Q recurrence (ambient n-1=15), specialized at (x,x,1).
    g = fmpq_poly([1,0,-1]); z = g*g
    qs = [fmpq_poly([1]), g]
    for k in range(1,16):
        qs.append(((2*k+13)*g*qs[-1] - k*z*qs[-2])/(k+13))
    assert all(q == g**k for k,q in enumerate(qs))
    pair = f + 2*B[0][1] + B[1][1] + R(1)
    cross0 = None
    for k in range(17):
        L, lam, X = blocks[17+k]
        p = polys(16+2*k,16-k)
        factor_q = factor_quadratic(L,lam,p)
        dense_q = dense_quadratic(X,p)
        assert factor_q == dense_q
        pair += theta[k]*qs[k]*factor_q
        if k == 0:
            cross0 = factor_cross(L,lam,p)
            dense_cross = sum((sum(X[i][j] for i in range(17))*p[j] for j in range(17)), fmpq_poly([]))
            assert cross0 == dense_cross
            pair += 2*theta[0]*cross0
    PAIR = pair
    C = B[0][0] + theta[0]*cross0(1)
    objective = f(1) + C
    source_objective = Q(1)+sum(blocks[k][2][0][0] for k in range(16))+B[0][0]+cross0(1)
    assert source_objective > 0 # New row weights; old objective not reused.
    assert Q(0) < objective < Q(8100)
    pair_payload = {
        'coefficient_order': 'ascending monomial powers x^0 through x^degree',
        'degree': pair.degree(),
        'coefficients': [str(v) for v in pair.coeffs()],
        'f_coefficients': [str(v) for v in f.coeffs()],
        'R_coefficients': [str(v) for v in R.coeffs()],
        'origin_C_exact': str(C),
        'objective_exact': str(objective),
    }
    dump('pair_coefficients.json', pair_payload)
    print('Exact reconstruction complete:', pair.degree(), decimal_enclosure(objective), flush=True)
    cert0 = sturm_certificate(pair,left,cap,'pair')
    print('Strict pair negativity proved.', flush=True)
    # Quantitative safe slack well below the approximately 9.99e-7 sampled slack.
    margin = Q('9/10000000')
    certm = sturm_certificate(pair+margin,left,cap,'pair_plus_9e_minus_7')
    result = {
        'status': 'PASS: exact PSD, nonnegative coefficients, origin/objective, and entire pair interval',
        'not_a_count_or_density_certificate': True,
        'missing_gate': 'No cyclic three-point inequality is checked here. The separate cyclic gate is indispensable.',
        'input_sha256': EXPECTED,
        'arithmetic': 'python-flint fmpq and fmpq_poly, exact rational decimals, threads=1',
        'dimension': 16,
        'cap_exact': str(cap),
        'interval_exact': [str(left),str(cap)],
        'source_factor_records': records,
        'all_source_lambdas_positive': True,
        'psd_proof': 'For each block X=L^T L+lambda I, z^T X z=||Lz||^2+lambda||z||^2 >= 0. All 34 lambda values are exact positive rationals. B is unchanged; theta_k F_k is PSD for every nonnegative exact theta_k.',
        'all_theta_a_r_nonnegative': True,
        'positive_local_weights': [k for k,v in enumerate(theta) if v > 0],
        'positive_global_coefficients': [j+1 for j,v in enumerate(a) if v > 0],
        'positive_nonorigin_R_coefficients': [j+1 for j,v in enumerate(r) if v > 0],
        'R_is_used_only_for_nonorigin_pairs': True,
        'source_objective_exactly_recomputed': True,
        'factor_and_dense_matrix_pair_reconstructions_agree': True,
        'homogeneous_Q_pair_recurrence_equals_one_minus_x_squared_power': True,
        'f1_exact': str(f(1)),
        'R1_exact': str(R(1)),
        'origin_C_exact': str(C),
        'origin_C_enclosure': decimal_enclosure(C),
        'objective_exact': str(objective),
        'objective_enclosure': decimal_enclosure(objective),
        'pair_polynomial_degree': pair.degree(),
        'pair_coefficients_sha256': sha((HERE/'pair_coefficients.json').read_bytes()),
        'pair_sturm': cert0,
        'strict_uniform_pair_upper_bound_exact': str(-margin),
        'pair_with_margin_sturm': certm,
        'runtime_seconds': monotonic()-start,
    }
    dump('result.json',result)
    print(json.dumps({k:result[k] for k in ['status','objective_enclosure','pair_polynomial_degree','strict_uniform_pair_upper_bound_exact','runtime_seconds']}, indent=2), flush=True)

if __name__ == '__main__':
    main()
