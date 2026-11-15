"""Exact periodic 93--96D existence bounds from the classical MW92 lattice.

Uses the repo's existing greedy-torus/colored-layer method. No explicit
horizontal offsets, world-priority claim, or new construction principle.
"""
from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import product
from math import factorial
from pathlib import Path
import argparse
import json
import sys

if not __debug__:
    raise RuntimeError("Run without -O.")
P = Path(__file__).resolve().parent
sys.path.insert(0, str(P.parent/"cross_dimension_lower_20261115"))
from verify_counting import root_gram, determinant, decsqrt
from verify_to96 import pi_interval


def mul(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def transpose(a):
    return list(map(list, zip(*a)))


def color_data(r):
    Q = root_gram(r)
    if r == 1:
        C, N, lam = [[2]], 2, F(4)
    elif r == 2:
        C, N, lam = [[1, -1], [1, 2]], 3, F(3)
    elif r == 3:
        # Kernel of a+b+2c mod 6. All nonzero vectors of norm below 6
        # are excluded exactly; the coefficient box is justified by Q^-1.
        C, N, lam = [[6, 0, 0], [-1, 1, 0], [-2, 0, 1]], 6, F(3)
        assert all(sum(a*b for a, b in zip(row, [1, 1, 2])) % 6 == 0 for row in C)
        Qinv = [[F(3, 4), F(1, 2), F(1, 4)], [F(1, 2), F(1), F(1, 2)],
                [F(1, 4), F(1, 2), F(3, 4)]]
        assert mul(Q, Qinv) == [[int(i == j) for j in range(3)] for i in range(3)]
        inverse_diagonal = [Qinv[i][i] for i in range(3)]
        assert all(6*x < 3**2 for x in inverse_diagonal)
        for z in product(range(-3, 4), repeat=3):
            norm = sum(z[i]*Q[i][j]*z[j] for i in range(3) for j in range(3))
            if 0 < norm < 6:
                assert (z[0]+z[1]+2*z[2]) % 6
        assert sum(C[1][i]*Q[i][j]*C[1][j] for i in range(3) for j in range(3)) == 6
    else:
        assert r == 4
        # Multiplication by 1+i+j in the quaternion model of D4 has norm 3.
        R = [[1, -1, 0, 0], [0, 1, -1, 0], [0, 0, 1, -1], [0, 0, 1, 1]]
        A = [[1, -1, -1, 0], [1, 1, 0, 1], [1, 0, 1, -1], [0, -1, 1, 1]]
        Rinv = [[F(1), F(1), F(1, 2), F(1, 2)],
                [F(0), F(1), F(1, 2), F(1, 2)],
                [F(0), F(0), F(1, 2), F(1, 2)],
                [F(0), F(0), F(-1, 2), F(1, 2)]]
        assert mul(R, Rinv) == [[int(i == j) for j in range(4)] for i in range(4)]
        C = mul(mul(R, A), Rinv)
        assert all(x.denominator == 1 for row in C for x in row)
        C = [[int(x) for x in row] for row in C]
        N, lam = 9, F(3)
    assert abs(determinant(C)) == N
    period = mul(mul(C, Q), transpose(C))
    if r != 3:
        assert period == [[lam*x for x in row] for row in Q]
    assert all(determinant([row[:j] for row in Q[:j]]) > 0 for j in range(1, r+1))
    return Q, C, N, lam


def audit():
    lo, pi_hi = pi_interval()
    assert lo < pi_hi
    parent_delta2 = F(2**184, 47**14)
    parent_volume = 47**7
    rows = []
    for r in range(1, 5):
        Q, C, N, lam = color_data(r)
        # Largest micro-unit depth satisfying the STRICT sufficient volume
        # inequality with a rigorous upper bound for pi.
        def fits(d):
            return (N-1)*pi_hi**46*d**46 < parent_volume*factorial(46)
        left, right, denominator = 0, 12_000_000, 1_000_000
        while left+1 < right:
            mid = (left+right)//2
            if fits(F(mid, denominator)):
                left = mid
            else:
                right = mid
        depth = F(left, denominator)
        assert fits(depth) and not fits(depth+F(1, denominator))
        vertical_minimum = 16-depth
        assert depth > 0 and vertical_minimum > 0
        assert lam*vertical_minimum >= 16
        vertical_det = (vertical_minimum/2)**r*determinant(Q)
        delta2 = parent_delta2*F(4)**r/vertical_det
        product2 = parent_delta2/F(2**r*determinant(Q))
        assert delta2/product2 == (16/vertical_minimum)**r
        # Same-color vectors use the period; different colors use horizontal
        # separation depth and the full vertical minimum.
        assert depth+vertical_minimum == 16
        row = {"dimension": 92+r, "colors": N,
               "vertical_root_gram": Q, "color_period_matrix": C,
               "same_color_squared_norm_ratio": str(lam),
               "horizontal_squared_distance_guarantee": str(depth),
               "vertical_squared_minimum": str(vertical_minimum),
               "strict_forbidden_volume_inequality_checked": True,
               "minimum_squared_exact_inheriting_MW92": 16,
               "period_cosets": N,
               "center_density_squared_exact": str(delta2),
               "center_density": decsqrt(delta2),
               "MW92_product_center_density_squared_exact": str(product2),
               "ratio_to_product_squared_exact": str(delta2/product2),
               "ratio_to_product": decsqrt(delta2/product2),
               "explicit_horizontal_offsets": False,
               "world_record_confirmed": False}
        rows.append(row)
    # BCS 1982, Eq.(15), m=48, n=2, a=1, parent delta=(3/2)^24.
    bcs96 = F(3**48, 2**24)
    assert bcs96 == F(3, 2)**48*2**24
    ratio2 = F(rows[-1]["center_density_squared_exact"])/bcs96**2
    with localcontext() as ctx:
        ctx.prec = 50
        ratio = (Decimal(ratio2.numerator)/Decimal(ratio2.denominator)).sqrt()
        gain = 100*(ratio-1)
        floor = int(gain*10**6)
    interval = [str(F(floor, 10**6)), str(F(floor+1, 10**6))]
    assert (1+F(interval[0])/100)**2 < ratio2 < (1+F(interval[1])/100)**2
    return {"dimension_ceiling": 96,
            "parent_input": "Shioda MW92: minimum 16, Gram determinant 47^14; Schuett-Shioda 2017 draft, Thm 13.35, p=47",
            "method_credit": "Classical greedy torus offsets and colored layers; already used in research/novel_routes_20261113/code_routes/TORUS_LAYER_BENCHMARKS.md",
            "rows": rows,
            "BCS96_comparison": {"source": "Bos-Conway-Sloane (1982), Theorem 2, Eq.(15), Table 1 and P48q example",
                                 "center_density_exact": "3^48/2^24", "center_density": decsqrt(bcs96**2),
                                 "center_density_squared_exact": str(bcs96**2),
                                 "new_to_BCS_ratio_squared_exact": str(ratio2),
                                 "new_to_BCS_ratio": str(ratio),
                                 "gain_percent": str(gain),
                                 "gain_percent_interval_exact": interval},
            "all_proof_inequalities_checked_exactly": True,
            "world_record_confirmed": False}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-result", action="store_true")
    args = ap.parse_args()
    result = audit()
    path = P/"torus93_96_results.json"
    if args.write_result:
        path.write_text(json.dumps(result, indent=2)+"\n")
    else:
        assert json.loads(path.read_text()) == result, "Stored result is stale."
    for row in result["rows"]:
        print(row["dimension"], row["colors"], row["center_density"], row["horizontal_squared_distance_guarantee"])
    print("96D gain over BCS1982:", result["BCS96_comparison"]["gain_percent"], "%")
