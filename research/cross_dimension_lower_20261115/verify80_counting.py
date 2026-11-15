"""Deeper compatible half-cosets above an extremal 80D parent.

Exact arithmetic for the extension in COUNTING_PROOF.md, Section 4.
The output proves existence, not explicit representatives or record status.
Standard library only. Reuses exact series and rational matrix operations.
"""
from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import argparse
import json
from verify_counting import theta, determinant, root_gram, decsqrt

if not __debug__:
    raise RuntimeError("Assertions must be enabled.")

P = Path(__file__).resolve().parent


def root_model(gram):
    r = len(gram)
    counts = [0, 0]
    # Graph in r hyperbolic planes: a -> (a,U*a). U has the quadratic
    # polynomial coefficients on/above its diagonal. This proves an embedding
    # for the ACTUAL restricted quadratic form, not merely matching counts.
    U = [[int(gram[i][i]//2) if i == j else gram[i][j] % 2 if i < j else 0
          for j in range(r)] for i in range(r)]
    for a in product((0, 1), repeat=r):
        norm = sum(a[i]*gram[i][j]*a[j] for i in range(r) for j in range(r))
        assert norm % 2 == 0
        q = (norm//2) % 2
        y = [sum(U[i][j]*a[j] for j in range(r)) % 2 for i in range(r)]
        assert sum(a[i]*y[i] for i in range(r)) % 2 == q
        if any(a):
            counts[q] += 1
    return counts


def audit():
    t, coeffs = theta(80)
    bad0 = sum(t[n//2] for n in (8, 12, 16))//2  # singular depth >=5
    bad1 = sum(t[n//2] for n in (10, 14, 18))//2  # nonsingular depth >=11/2
    population = 2**79-2**39-1
    rows = []
    for r in range(1, 6):
        # D5's actual quadratic form has 20 nonsingular vectors and FAILS
        # this union bound. A5 has 16 and proves the 85D case.
        root = root_gram(r, family="A" if r == 5 else None)
        z, o = root_model(root)
        assert [z, o] == [[0, 1], [0, 3], [1, 6], [3, 12], [15, 16]][r-1]
        expect = F(z*bad0+o*bad1, population)
        assert expect < 1
        gram = [[F(5, 4)*x for x in row] for row in root]
        det = determinant(gram)
        rootdet = [2, 3, 4, 4, 6][r-1]
        assert det == F(5, 4)**r*rootdet
        assert all(determinant([row[:j] for row in gram[:j]]) > 0
                   for j in range(1, r+1))
        # Root lattice norms are positive even integers. q=1 permits norm 2,
        # q=0 forces norm >=4; after scaling these are >=5/2 and >=5.
        assert F(11, 2)+F(5, 2) == 8
        assert 5+5 >= 8
        assert 4*F(5, 2) >= 8  # all-even coefficients, nonzero 2M vector
        d2 = F(2**(80+r), 1)/det
        comparison_det = [2, 3, 4, 4, 4][r-1]
        product2 = F(2**80, 2**r*comparison_det)
        ratio2 = d2/product2
        assert ratio2 == F(16, 5)**r*F(comparison_det, rootdet)
        rows.append({
            "dimension": 80+r, "rank": r,
            "nonzero_singular_classes_in_model": z,
            "nonsingular_classes_in_model": o,
            "expected_bad_classes_upper_bound_exact": str(expect),
            "vertical_gram_exact": [[str(x) for x in row] for row in gram],
            "vertical_root_lattice": ["A1", "A2", "A3", "D4", "A5"][r-1],
            "gram_determinant_exact": str(det),
            "minimum_squared_exact_inheriting_parent_minimum": 8,
            "center_density_squared_exact": str(d2),
            "center_density": decsqrt(d2),
            "comparison_method": "Bachoc--Nebe L80 orthogonal product with minimum-2 A1/A2/A3/D4/D5, rescaled to minimum 8",
            "comparison_is_reported_dimension_specific_record": False,
            "product_center_density_squared_exact": str(product2),
            "ratio_to_product_squared_exact": str(ratio2),
            "ratio_to_product": decsqrt(ratio2),
            "explicit_horizontal_coset_representatives": False,
            "published_record_improvement_established": False,
        })
    return {
        "status": "exact_arithmetic_checked_for_80D_parent_existence_proof",
        "parent_input": "Bachoc--Nebe extremal even unimodular L80 of minimum 8; the same argument works for M80",
        "source_url": "https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/L80.html",
        "proof": "COUNTING_PROOF.md#4-the-80d-parent-gives-deeper-nonsingular-half-cosets",
        "theta_modular_basis_coefficients": coeffs,
        "shell_counts_by_squared_norm": {str(2*j): t[j] for j in range(4, len(t))},
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
    path = P/"counting80_results.json"
    if args.write_result:
        path.write_text(json.dumps(result, indent=2)+"\n")
    else:
        assert json.loads(path.read_text()) == result, "Stored 80D-parent report is stale."
    for row in result["rows"]:
        print(f'{row["dimension"]}D: determinant={row["gram_determinant_exact"]}, '
              f'center density={row["center_density"]}, '
              f'ratio to named product={row["ratio_to_product"]}')
