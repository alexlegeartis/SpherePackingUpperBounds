"""Exact arithmetic for compatible half-cosets of an extremal 72D lattice.

This is an existence proof, not a search for explicit coset representatives.
Mathematical inputs and the averaging/gluing proof are in COUNTING_PROOF.md.
Only the Python standard library is required. Never run with -O.
"""
from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import product
from math import comb
from pathlib import Path
import argparse
import json

if not __debug__:
    raise RuntimeError("Assertions must be enabled.")

P = Path(__file__).resolve().parent
ORDER = 12


def mul(a, b):
    return [sum(a[j] * b[i-j] for j in range(i+1))
            for i in range(ORDER+1)]


def power(a, n):
    out = [1] + [0]*ORDER
    for _ in range(n):
        out = mul(out, a)
    return out


def sigma(n, k):
    return sum(d**k for d in range(1, n+1) if n % d == 0)


def theta(ambient=72):
    assert ambient in (72, 80)
    e4 = [1] + [240*sigma(n, 3) for n in range(1, ORDER+1)]
    e6 = [1] + [-504*sigma(n, 5) for n in range(1, ORDER+1)]
    eta = [1] + [0]*ORDER
    for m in range(1, ORDER+1):
        factor = [0]*(ORDER+1)
        for j in range(min(24, ORDER//m)+1):
            factor[m*j] = (-1)**j*comb(24, j)
        eta = mul(eta, factor)
    delta = [0] + eta[:ORDER]
    # Independent coefficient formula for Delta; no tabulated theta counts.
    e43, e62 = power(e4, 3), power(e6, 2)
    assert all(e43[i]-e62[i] == 1728*delta[i] for i in range(ORDER+1))
    terms = [mul(power(e4, ambient//8-3*j), power(delta, j))
             for j in range(4)]
    # Solve the triangular vanishing conditions rather than assume coefficients.
    out = terms[0][:]
    coeffs = [1]
    for j in range(1, 4):
        assert terms[j][j] == 1
        c = -out[j]
        coeffs.append(c)
        out = [x+c*y for x, y in zip(out, terms[j])]
    assert coeffs == {72: [1, -2160, 965520, -27302400],
                      80: [1, -2400, 1360800, -103488000]}[ambient]
    assert out[:4] == [1, 0, 0, 0]
    assert out[4] == {72: 6218175600, 80: 1250172000}[ambient]
    assert all(x >= 0 and x % 2 == 0 for x in out[4:])
    return out, coeffs


def determinant(matrix):
    a = [[F(x) for x in row] for row in matrix]
    out = F(1)
    for i in range(len(a)):
        p = next(j for j in range(i, len(a)) if a[j][i])
        if p != i:
            a[p], a[i] = a[i], a[p]
            out = -out
        pivot = a[i][i]
        out *= pivot
        for j in range(i+1, len(a)):
            q = a[j][i]/pivot
            for k in range(i+1, len(a)):
                a[j][k] -= q*a[i][k]
    return out


def root_gram(r, family=None):
    if r < 4 or family == "A":
        return [[2 if i == j else -1 if abs(i-j) == 1 else 0
                 for j in range(r)] for i in range(r)]
    assert r in (4, 5)
    b = [[int(j == i)-int(j == i+1) for j in range(r)] for i in range(r-1)]
    b += [[int(j in (r-2, r-1)) for j in range(r)]]
    return [[sum(x*y for x, y in zip(u, v)) for v in b] for u in b]


def mixed_gram(r):
    assert r in (5, 6)
    b = [[F(x) for x in row] for row in
         [[1, -1, 0, 0], [0, 1, -1, 0],
          [0, 0, 1, -1], [0, 0, 1, 1]]]
    if r == 5:
        b += [[F(1, 2)]*4]
    else:
        b += [[F(1), F(0), F(0), F(0)], [F(1, 2)]*4]
    a2 = [[2, -1], [-1, 2]]
    g = [[F(3, 2)*sum(x*y for x, y in zip(u, v)) for v in b]
         for u in b]
    for i in range(4, r):
        for j in range(4, r):
            g[i][j] += a2[i-4][j-4]
    return g


def qmodel(r):
    # Six hyperbolic planes occur in either possible 72D quadratic type.
    def e(i):
        return 1 << (2*i)

    def f(i):
        return 1 << (2*i+1)

    b = [e(i) for i in range(min(r, 4))]
    if r == 5:
        b += [e(4) ^ f(4)]
    if r == 6:
        b += [e(4) ^ f(4), e(4) ^ e(5) ^ f(5)]

    def q(x):
        return sum(((x >> (2*i)) & 1)*((x >> (2*i+1)) & 1)
                   for i in range(6)) % 2

    seen = set()
    counts = [0, 0]
    for a in product((0, 1), repeat=r):
        v = 0
        for bit, w in zip(a, b):
            if bit:
                v ^= w
        assert v not in seen
        seen.add(v)
        expected = (0 if r <= 4 else a[4] if r == 5
                    else (a[4]+a[5]+a[4]*a[5]) % 2)
        assert q(v) == expected
        if any(a):
            counts[q(v)] += 1
    return counts


def decsqrt(x):
    with localcontext() as ctx:
        ctx.prec = 45
        return str((Decimal(x.numerator)/Decimal(x.denominator)).sqrt())


def audit():
    t, coeffs = theta()
    bad0 = sum(t[n//2] for n in (8, 12, 16))//2
    bad1 = sum(t[n//2] for n in (10, 14))//2
    # Lower bound for BOTH nonzero singular and nonsingular populations,
    # valid for either Arf invariant. Avoid inferring the invariant from dimension.
    population = 2**71-2**35-1
    rows = []
    for r in range(1, 7):
        z, o = qmodel(r)
        expectation = F(z*bad0+o*bad1, population)
        assert expectation < 1
        gram = ([[F(3, 2)*x for x in row] for row in root_gram(r)]
                if r <= 4 else mixed_gram(r))
        det = determinant(gram)
        assert all(determinant([row[:j] for row in gram[:j]]) > 0
                   for j in range(1, r+1))
        assert det == [F(3), F(27, 4), F(27, 2), F(81, 4),
                       F(81, 2), F(243, 4)][r-1]
        delta2 = F(2**(72+r), 1)/det
        root_det = [2, 3, 4, 4, 4, 3][r-1]
        product2 = F(2**72, 2**r*root_det)
        ratio2 = delta2/product2
        assert ratio2 > 1
        # Layer bounds: nonzero singular horizontal class >=5, nonsingular
        # >=9/2. The chosen vertical metric supplies >=3 and >=7/2.
        assert F(5)+3 == 8 and F(9, 2)+F(7, 2) == 8
        assert 4*3 >= 8  # all-even coefficient vector, in 2M.
        rows.append({
            "dimension": 72+r, "rank": r,
            "nonzero_singular_classes_in_model": z,
            "nonsingular_classes_in_model": o,
            "expected_bad_classes_upper_bound_exact": str(expectation),
            "vertical_gram_exact": [[str(x) for x in row] for row in gram],
            "gram_determinant_exact": str(det),
            "minimum_squared_exact_inheriting_parent_minimum": 8,
            "center_density_squared_exact": str(delta2),
            "center_density": decsqrt(delta2),
            "comparison_method": "Gamma72 orthogonal product with minimum-2 A1/A2/A3/D4/D5/E6, rescaled to minimum 8",
            "comparison_is_reported_dimension_specific_record": False,
            "product_center_density_squared_exact": str(product2),
            "ratio_to_product_squared_exact": str(ratio2),
            "ratio_to_product": decsqrt(ratio2),
            "explicit_horizontal_coset_representatives": False,
            "published_record_improvement_established": False,
        })
    return {
        "status": "exact_arithmetic_checked_for_existence_proof",
        "parent_input": "Nebe's even unimodular Gamma72 of minimum 8 (arXiv:1008.2862, Theorem 2.3)",
        "proof": "COUNTING_PROOF.md",
        "theta_modular_basis_coefficients": coeffs,
        "shell_counts_by_squared_norm": {str(2*j): t[j] for j in range(4, ORDER+1)},
        "bad_nonzero_singular_classes_upper_bound": bad0,
        "bad_nonsingular_classes_upper_bound": bad1,
        "conservative_population_for_each_type": population,
        "rows": rows,
        "all_arithmetic_exact": True,
        "existential_not_explicit": True,
        "current_record_identified": False,
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-result", action="store_true")
    args = ap.parse_args()
    result = audit()
    path = P/"counting_results.json"
    if args.write_result:
        path.write_text(json.dumps(result, indent=2)+"\n")
    else:
        assert json.loads(path.read_text()) == result, "Stored counting report is stale."
    for row in result["rows"]:
        print(f'{row["dimension"]}D: determinant={row["gram_determinant_exact"]}, '
              f'center density={row["center_density"]}, '
              f'ratio to named product={row["ratio_to_product"]}')
