#!/usr/bin/env python3
"""Independent exact Ne64 theta and discriminant-count verification."""

from pathlib import Path
from math import comb
import json

import numpy as np


HERE = Path(__file__).resolve().parent
DEGREE = 13


def mul(a, b):
    c = [0] * (DEGREE + 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b[:DEGREE + 1 - i]):
            c[i + j] += x * y
    return c


def power(a, n):
    c = [1] + [0] * DEGREE
    for _ in range(n):
        c = mul(c, a)
    return c


def a2_theta():
    # Divisor formula, independent of direct enumeration of A2 vectors.
    a = [1]
    for n in range(1, DEGREE + 1):
        a.append(6 * sum(0 if d % 3 == 0 else (1 if d % 3 == 1 else -1)
                         for d in range(1, n + 1) if n % d == 0))
    return a


def eta_product():
    a = [0, 1] + [0] * (DEGREE - 1)
    for n in range(1, DEGREE + 1):
        for exponent in (n, 3 * n):
            b = [0] * (DEGREE + 1)
            for j in range(7):
                if exponent * j <= DEGREE:
                    b[exponent * j] = (-1)**j * comb(6, j)
            a = mul(a, b)
    return a


def row_basis_mod3(a):
    rows = []
    basis = {}
    for i, row in enumerate(a):
        v = row.copy() % 3
        for pivot, b in sorted(basis.items()):
            v = (v - v[pivot] * b) % 3
        if np.any(v):
            pivot = int(np.flatnonzero(v)[0])
            v = v * int(v[pivot]) % 3
            basis[pivot] = v
            rows.append(i)
    return rows


def determinant_mod3(a):
    a = a.copy() % 3
    det = 1
    for j in range(len(a)):
        p = next((i for i in range(j, len(a)) if a[i, j]), None)
        if p is None:
            return 0
        if p != j:
            a[[p, j]] = a[[j, p]]
            det = -det
        pivot = int(a[j, j])
        det *= pivot
        a[j] = a[j] * pivot % 3
        for i in range(j + 1, len(a)):
            a[i] = (a[i] - a[i, j] * a[j]) % 3
    return det % 3


def main():
    a, d = a2_theta(), eta_product()
    forms = [mul(power(a, 32 - 6*j), power(d, j)) for j in range(6)]
    series = [0] * (DEGREE + 1)
    c = []
    for j in range(6):
        assert forms[j][:j] == [0] * j and forms[j][j] == 1
        coefficient = (1 if j == 0 else 0) - series[j]
        c.append(coefficient)
        series = [x + coefficient * y for x, y in zip(series, forms[j])]
    assert c == [1, -192, 10944, -185088, 536256, -483840]
    assert series[1:6] == [0] * 5
    assert series[6] == 138458880
    assert series[7] == 15482327040
    assert series[9] == 37527682805760
    assert series[10] == 983531592449280

    H = np.load(HERE.parent / "dual64" / "source.npz")["H"]
    rows = row_basis_mod3(H)
    assert len(rows) == 32
    # Q(c)=3 (c,c)/2 mod3, and c_i=e_i G^{-1}; its matrix is H/2.
    gram = 2 * H[np.ix_(rows, rows)] % 3
    determinant = determinant_mod3(gram)
    assert determinant == 1
    # In even dimension 2m, plus type iff (-1)^m det is a square.
    assert len(rows) // 2 % 2 == 0
    all_isotropic = 3**31 + 2 * 3**15
    nonzero_isotropic = all_isotropic - 1
    each_nonisotropic = 3**31 - 3**15
    bad_isotropic_upper_bound = series[6] + series[9]
    good_isotropic_lower_bound = nonzero_isotropic - bad_isotropic_upper_bound
    assert good_isotropic_lower_bound > 0
    norm14_and_20_count = series[7] + series[10]
    assert norm14_and_20_count > each_nonisotropic
    assert norm14_and_20_count // 3 < each_nonisotropic

    result = {
        "theta_basis_coefficients": c,
        "theta_coefficients_q_power_equals_half_norm": series,
        "discriminant_basis_indices": rows,
        "discriminant_quadratic_matrix_determinant_mod3": determinant,
        "discriminant_type": "plus",
        "nonzero_isotropic_classes": nonzero_isotropic,
        "each_nonisotropic_value_classes": each_nonisotropic,
        "isotropic_classes_representable_below_norm8_upper_bound": bad_isotropic_upper_bound,
        "isotropic_classes_of_minimum_at_least8_lower_bound": good_isotropic_lower_bound,
        "nonisotropic_short_vector_count_norm14_over3_or20_over3": norm14_and_20_count,
        "naive_count_proves_norm26_over3_hole": False,
        "order_three_multiplicity_required_but_not_established": True,
    }
    (HERE / "ne64_theta_count_verified.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
