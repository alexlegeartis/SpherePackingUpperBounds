"""Reconstruct an explicit six-translate packing from its primitive normal.

This writes exact integer/rational data. The separate verifier re-parses the
published ambient and checks all identities without importing this program.
"""
from pathlib import Path
import itertools
import math
import sympy as sp

P = Path(__file__).resolve().parent
ROOT = P.parents[1]


def read(name):
    return sp.Matrix([[int(x) for x in row.split()]
                      for row in (P / "data" / name).read_text().splitlines()])


def write(name, M):
    assert all(x.q == 1 for x in M)
    (P / "data" / name).write_text("\n".join(" ".join(map(str, M.row(i)))
                                                for i in range(M.rows)) + "\n")


def main():
    # The archived rational dual identity is used only by this data generator.
    # verify.py obtains G and its inverse independently from catalogue HTML.
    import numpy as np
    source = np.load(ROOT / "research/novel_routes_20261113/dual64/source.npz")
    G = sp.Matrix(source["G"].tolist())
    Gi = G.inv()
    A, S = read("normal_coefficients.txt"), read("dual_cluster.txt")
    K = A * Gi * A.T
    assert 3*K == read("normal_gram_numerator.txt")
    for cols in itertools.combinations(range(64), 5):
        C = A[:, list(cols)]
        if abs(C.det()) == 1:
            unit = list(cols)
            break
    free = [i for i in range(64) if i not in unit]
    Ci = C.inv()
    right = sp.zeros(64, 5)
    for j, c in enumerate(unit): right[c, :] = Ci[j, :]
    N = sp.zeros(59, 64)
    left = -A[:, free].T * Ci.T
    for j, c in enumerate(unit): N[:, c] = left[:, j]
    for j, c in enumerate(free): N[j, c] = 1
    lifts = S * right.T
    projected = lifts - S * K.inv() * A * Gi
    shifts = projected[:, free]
    den = math.lcm(*(int(v.q) for v in shifts))
    write("unit_minor_columns.txt", sp.Matrix([unit]))
    write("section_basis.txt", N)
    write("period_gram.txt", N * G * N.T)
    write("ambient_lifts.txt", lifts)
    write("shift_numerators.txt", shifts * den)
    write("shift_denominator.txt", sp.Matrix([[den]]))
    print("Wrote explicit 59D period and six translates; denominator", den)


if __name__ == "__main__":
    main()
