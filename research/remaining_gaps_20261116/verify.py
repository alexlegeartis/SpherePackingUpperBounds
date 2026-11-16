"""Exact certificates from the continuation on the nine unfilled dimensions.

This package claims no new sphere-packing record. Published minimum and
3-modularity of Ne64 are inherited inputs, not re-proved by enumeration.
"""
from pathlib import Path
from fractions import Fraction
import argparse
import hashlib
from itertools import combinations, permutations
import json
import re

import numpy as np
from flint import fmpz_mat, fmpq_mat, nmod_mat, fmpz_mod_poly_ctx, fmpz_poly

if not __debug__:
    raise RuntimeError("Run without -O.")

P = Path(__file__).resolve().parent
OLD = P.parent / "novel_routes_20261113"


def source_gram(path):
    raw = path.read_text()
    section = raw.split('<a NAME="GRAM">')[1].split('<a NAME="DIVISORS">')[0]
    values = list(map(int, re.findall(r"-?\d+", re.sub("<[^>]*>", " ", section))))
    assert values[:2] == [64, 0] and len(values) == 2 + 64 * 65 // 2
    rows = [[0] * 64 for _ in range(64)]
    pos = 2
    for i in range(64):
        for j in range(i + 1):
            rows[i][j] = rows[j][i] = values[pos]
            pos += 1
    return fmpz_mat(rows)


def integral(Q):
    assert all(x.denominator == 1 for row in Q.tolist() for x in row)
    return fmpz_mat([[int(x) for x in row] for row in Q.tolist()])


def ne64_geometry():
    path = OLD / "dual64/CQ64b.html"
    raw = path.read_text()
    minimum_section = raw.split('<a NAME="MINIMAL_NORM">')[1].split("<p><li>")[0]
    assert re.sub("<[^>]*>", " ", minimum_section).split() == ["MINIMAL_NORM", "12"]
    assert "MODULAR=3" in raw
    G = source_gram(path)
    assert G == G.transpose() and G.det() == 3**32
    assert all(G[i, i] % 2 == 0 for i in range(64))
    H = integral(3 * G.inv())
    assert all(H[i, i] % 2 == 0 for i in range(64))
    # Positive definiteness, without a numerical eigensolver.
    assert all(fmpz_mat([r[:k] for r in G.tolist()[:k]]).det() > 0
               for k in range(1, 65))
    with np.load(P / "ne64_perfection.npz", allow_pickle=False) as data:
        V = data["V"].astype(np.int64)
    assert V.shape == (2080, 64)
    A = np.array(G.tolist(), dtype=np.int64)
    assert 64**2 * int(abs(V).max())**2 * int(abs(A).max()) < 2**63
    assert np.all(np.einsum("vi,ij,vj->v", V, A, V) == 12)
    i, j = np.triu_indices(64)
    rank = nmod_mat((V[:, i] * V[:, j]).tolist(), 65521).rank()
    assert rank == 64 * 65 // 2
    # The harmonic-theta argument is given in PROOF.md. These are the
    # exact valence arithmetic and hypotheses used there.
    harmonic_weight = 64 // 2 + 2
    cusp_order_sum = 6 + 6
    valence_bound = Fraction(harmonic_weight * 4, 12)
    assert cusp_order_sum > valence_bound
    return {
        "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "inherited_published_inputs": {"minimum": 12, "modularity": 3},
        "certificate_rays": len(V), "dyad_rank_mod65521": rank,
        "harmonic_theta_weight": harmonic_weight,
        "harmonic_theta_cusp_order_sum_lower_bound": cusp_order_sum,
        "harmonic_theta_valence_upper_bound": str(valence_bound),
        "perfection_verified_exactly": True,
        "local_lattice_optimum_with_published_inputs_and_theorems": True,
        "global_optimality_claimed": False,
    }


def form_family():
    G = source_gram(OLD / "dual64/CQ64b.html")
    G0 = source_gram(OLD / "code_routes/CQ64a.html")
    assert G0.det() == 1
    assert all(fmpz_mat([r[:k] for r in G0.tolist()[:k]]).det() > 0
               for k in range(1, 65))
    Gi = integral(G0.inv())
    forms = [G0, G, G * Gi * G, G * Gi * G * Gi * G]
    ij = [(i, j) for i in range(64) for j in range(i, 64)]
    assert all(g[i, i] % 2 == 0 for g in forms for i in range(64))
    vectors = [[int(g[i, j] // 2) if i == j else int(g[i, j])
                for i, j in ij] for g in forms]
    assert fmpz_mat(vectors).rank() == 4
    # HNF of the coordinate functionals determines the coefficient lattice
    # giving integral forms with even diagonal, throughout this Q-span.
    A = fmpz_mat(vectors).transpose().hnf()
    assert all(x == 0 for row in A.tolist()[4:] for x in row)
    K = fmpz_mat(A.tolist()[:4])
    I = K.inv()
    base = []
    for i in range(4):
        f = fmpq_mat(64, 64)
        for j, g in enumerate(forms):
            f += I[j, i] * fmpq_mat(g)
        base.append(integral(f))
    with np.load(P / "global64_forms.npz", allow_pickle=False) as data:
        saved = data["forms"]
    assert saved.shape == (4, 64, 64)
    assert [fmpz_mat(a.tolist()) for a in saved] == base
    assert all(f[i, i] % 2 == 0 for f in base for i in range(64))
    with np.load(OLD / "code_routes/cq64_generators.npz", allow_pickle=False) as data:
        maps = [fmpz_mat(a.tolist()) for a in data["generators"][:5]]
    assert all(abs(C.det()) == 1 and C * f * C.transpose() == f
               for C in maps for f in base)
    verify_norm_polynomials(G0, base)
    return {"rational_family_dimension": 4,
            "unimodular_source_sha256": hashlib.sha256(
                (OLD / "code_routes/CQ64a.html").read_bytes()).hexdigest(),
            "even_integral_saturation_hnf": [[int(x) for x in r] for r in K.tolist()],
            "integral_basis_reconstructed_exactly": True,
            "entire_invariant_space_classification_claimed": False}


def polynomial_add(a, b, scale=1):
    out = a.copy()
    for key, value in b.items():
        out[key] = out.get(key, Fraction(0)) + scale * value
        if out[key] == 0:
            del out[key]
    return out


def polynomial_product(a, b):
    out = {}
    for u, x in a.items():
        for v, y in b.items():
            key = tuple(i + j for i, j in zip(u, v))
            out[key] = out.get(key, Fraction(0)) + x * y
    return {key: value for key, value in out.items() if value}


def verify_norm_polynomials(G0, base):
    """Rebuild the quartic field algebra and its characteristic coefficients."""
    data = json.loads((P / "norm_polynomials.json").read_text())
    T = [G0.inv() * fmpq_mat(f) for f in base]
    assert T[0] == fmpq_mat(fmpz_mat(np.eye(64, dtype=int).tolist()))
    pos = data["positions"]
    A = fmpq_mat([[a[i, j] for a in T] for i, j in pos])
    Ai = A.inv()
    multiplication = []
    for a in T:
        columns = []
        for b in T:
            product = a * b
            c = Ai * fmpq_mat([[product[i, j]] for i, j in pos])
            rebuilt = fmpq_mat(64, 64)
            for k in range(4):
                rebuilt += c[k, 0] * T[k]
            assert rebuilt == product
            columns.append([c[k, 0] for k in range(4)])
        multiplication.append(fmpq_mat(columns).transpose())
    # The Q-algebra is a degree-four field; the 64D representation is
    # sixteen copies of its regular representation.
    cp = multiplication[1].charpoly()
    assert all(c.denominator == 1 for c in cp.coeffs())
    p = fmpz_poly([int(c) for c in cp.coeffs()])
    _, factors = p.factor()
    assert len(factors) == 1 and factors[0] == (p, 1) and p.degree() == 4
    assert integral(T[1]).charpoly() == p**16
    linear = [[{tuple(int(h == k) for h in range(4)):
                 Fraction(str(multiplication[k][i, j]))
                 for k in range(4) if multiplication[k][i, j]}
                for j in range(4)] for i in range(4)]
    zero = (0, 0, 0, 0)
    elementary = []
    for size in range(1, 5):
        total = {}
        for subset in combinations(range(4), size):
            for perm in permutations(range(size)):
                sign = (-1)**sum(perm[i] > perm[j]
                                 for i in range(size) for j in range(i + 1, size))
                term = {zero: Fraction(sign)}
                for i in range(size):
                    term = polynomial_product(term, linear[subset[i]][subset[perm[i]]])
                total = polynomial_add(total, term)
        elementary.append(total)
    saved = [{tuple(a["powers"]): Fraction(a["coefficient"]) for a in terms}
             for terms in data["elementary_symmetric_polynomials"]]
    assert elementary == saved


def code93_shift_obstruction():
    # x^31-1 is square-free in characteristic 2, so every invariant
    # subspace is a direct sum of simple F2[C31]-modules.
    R = fmpz_mod_poly_ctx(2)
    x = R([0, 1])
    f = x**31 - 1
    unit, factors = f.factor()
    assert unit == 1 and all(exponent == 1 for _, exponent in factors)
    degrees = sorted(g.degree() for g, _ in factors)
    assert degrees == [1, 5, 5, 5, 5, 5, 5]
    assert next(j for j in range(1, 31) if pow(2, j, 31) == 1) == 5
    # Three 31-cycles have a three-dimensional fixed subspace. Since 31
    # is odd, evenness cuts its dimension to at most two. Other simple
    # modules have dimension five. Therefore k mod5 is 0,1,or2.
    allowed = [0, 1, 2]
    assert 43 % 5 not in allowed
    return {"binary_length": 93, "dimension": 43,
            "irreducible_factor_degrees_x31_minus1": degrees,
            "even_fixed_subspace_dimension_upper": 2,
            "allowed_dimension_residues_mod5": allowed,
            "shift_by_three_invariant_even_code_impossible": True,
            "ordinary_cyclic_even_code_also_impossible": True,
            "general_even_binary_93_43_24_code_excluded": False}


def audit():
    from scan import scan
    return {"status": "no_new_packing_record", "new_packing_SOTA_claimed": False,
            "remaining_dimensions": [48, 54, 55, 56, 63, 64, 72, 80, 92],
            "ne64_geometry": ne64_geometry(), "ne64_form_family": form_family(),
            "bounded_ne64_family_search": scan(),
            "code93_shift_obstruction": code93_shift_obstruction()}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-result", action="store_true")
    args = ap.parse_args()
    result = audit()
    if args.write_result:
        (P / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    else:
        assert result == json.loads((P / "results.json").read_text())
    print("64D: perfection rank2080, 2237409 primitive forms rejected exactly; "
          "even 3-quasicyclic [93,43] impossible; no new packing SOTA")
