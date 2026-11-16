"""Exact finite rejection scan; no floating arithmetic or minimum oracle.

Run verify.py first to reconstruct the form basis and field polynomials.
The finite vector list is checked directly; its discovery is not trusted.
"""
from pathlib import Path
from math import isqrt
import json
import numpy as np

P = Path(__file__).resolve().parent


def scan(radius=32):
    assert isinstance(radius, int) and 1 <= radius <= 32
    with np.load(P / "global64_forms.npz", allow_pickle=False) as d:
        forms = d["forms"]
    with np.load(P / "scan64_witnesses.npz", allow_pickle=False) as d:
        V = d["V"].astype(np.int64)
    assert V.shape == (630, 64) and np.all(np.any(V != 0, axis=1))
    assert 64**2 * int(abs(V).max())**2 * int(abs(forms).max()) < 2**63
    Q = np.array([np.einsum("ij,ij->i", V @ f, V) for f in forms]).T
    assert 4 * radius * int(abs(Q).max()) < 2**63
    polynomials = json.loads((P / "norm_polynomials.json").read_text())[
        "elementary_symmetric_polynomials"]
    for terms in polynomials:
        assert sum(abs(a["coefficient"]) * radius**sum(a["powers"])
                   for a in terms) < 2**63
    count = positive = primitive = 0
    survivors = []
    # All tuples are visited. No approximate eigenvalue filter is used.
    tail = np.array([(b, c, d) for b in range(-radius, radius + 1)
                     for c in range(-radius, radius + 1)
                     for d in range(-radius, radius + 1)], dtype=np.int64)
    for a0 in range(-radius, radius + 1):
        Z = np.column_stack([np.full(len(tail), a0, dtype=np.int64), tail])
        count += len(Z)
        for terms in polynomials:
            value = np.zeros(len(Z), dtype=np.int64)
            for term in terms:
                values = np.full(len(Z), term["coefficient"], dtype=np.int64)
                for i, power in enumerate(term["powers"]):
                    values *= Z[:, i]**power
                value += values
            ok = value > 0
            Z, value = Z[ok], value[ok]
        # All four real eigenvalues are positive iff all four elementary
        # symmetric coefficients are positive. The last value is Norm(beta).
        norm = value
        positive += len(Z)
        ok = np.gcd.reduce(abs(Z), axis=1) == 1
        Z, norm = Z[ok], norm[ok]
        primitive += len(Z)
        threshold = np.array([isqrt(isqrt(2304 * int(n))) for n in norm],
                             dtype=np.int64)
        for q in Q:
            mu = Z @ q
            assert np.all(mu > 0)
            # delta(F)^2 <= (mu/4)^64/Norm(beta)^16.
            # Improvement over Ne64 requires mu^4 > 2304 Norm(beta).
            ok = mu > threshold
            Z, threshold = Z[ok], threshold[ok]
            if not len(Z):
                break
        survivors.extend(Z.tolist())
    assert count == (2 * radius + 1)**4
    if radius == 32:
        assert (positive, primitive) == (2427638, 2237409)
        assert survivors == []
    return {"coefficient_box_radius": radius,
            "integer_coefficient_tuples": count,
            "positive_definite_forms": positive,
            "primitive_positive_definite_forms": primitive,
            "integer_test_vectors": len(V), "not_rejected": survivors,
            "all_decisions_exact": True,
            "scope": "Only the supplied four-dimensional integral form family in this coefficient box; not a global density bound."}
