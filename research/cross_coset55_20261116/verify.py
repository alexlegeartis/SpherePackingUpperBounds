"""Exact cross-coset tilt in 55D, with finite negative certificates in 54-56D."""

from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext
import argparse
import json
import sys

import numpy as np
from flint import fmpz_mat, fmpq_mat, nmod_mat

if not __debug__:
    raise RuntimeError("Run without -O.")
P = Path(__file__).resolve().parent
sys.path.insert(0, str(P.parent / "cm58_20261116"))
import verify_frontier as frontier


def decsqrt(x):
    with localcontext() as ctx:
        ctx.prec = 70
        return str((Decimal(x.numerator) / x.denominator).sqrt())


def gain(ratio_squared):
    with localcontext() as ctx:
        ctx.prec = 70
        value = 100 * (Decimal(decsqrt(ratio_squared)) - 1)
        lo = F(int(value * 10**6), 10**6)
    hi = lo + F(1, 10**6)
    assert (1 + lo / 100) ** 2 <= ratio_squared < (1 + hi / 100) ** 2
    return {"percent": str(value), "interval": [str(lo), str(hi)]}


def build(G, cross, bottom):
    n = G.nrows()
    full = fmpq_mat(n + 1, n + 1)
    for i in range(n):
        for j in range(n):
            full[i, j] = G[i, j]
        full[i, n] = full[n, i] = cross[i]
    full[n, n] = bottom
    return full


def negative_audits(G54, q):
    bad = json.loads((P / "rejected_deformation55.json").read_text())
    gram = build(G54, bad["candidate_cross"], bad["candidate_bottom"])
    w = fmpq_mat([bad["exact_short_vector"]])
    norm = F(str((w * gram * w.transpose())[0, 0]))
    assert 0 < norm < 10 and str(norm) == bad["exact_squared_norm"]
    z = fmpq_mat([bad["new_algebraic_closest_vector_numerator"]]) / 3
    assert z * fmpq_mat(G54) * z.transpose() == fmpq_mat([["28/3"]])
    assert z - q == fmpq_mat([bad["exact_short_vector"][:54]])
    cross = fmpq_mat([bad["candidate_cross"]])
    h2 = F(bad["candidate_bottom"]) - F(
        str((cross * G54.inv() * cross.transpose())[0, 0])
    )
    assert h2 > 0
    # The failed proposal really had a smaller determinant, but not minimum 10.
    assert h2 < F(442, 405)
    assert F(442, 405) / h2 * (norm / 10) ** 55 < 1

    G56, _, _, _ = frontier.parent56()
    screen = json.loads((P / "overlattice_screen.json").read_text())
    assert screen["seed"] == 202611165456
    rows = screen["rows"]
    assert len(rows) == 48
    summaries = []
    for n, p, index, G in [(54, 3, 9, G54), (56, 29, 29, G56)]:
        trials = [a for a in rows if a["dimension"] == n]
        assert {a["trial"] for a in trials} == set(range(24))
        bounds = []
        for a in trials:
            assert a["prime"] == p and a["index"] == index
            r = fmpz_mat([a["added_numerator"]])
            assert all(int(x) % p == 0 for row in (r * G).tolist() for x in row)
            generators = (p * np.eye(n, dtype=np.int64)).tolist() + r.tolist()
            if n == 54:
                generators += [[int(x) for x in (3 * q).tolist()[0]]]
            B = fmpz_mat(fmpz_mat(generators).hnf().tolist()[:n])
            assert abs(B.det()) == p**n // index
            witness = fmpz_mat([a["witness_numerator"]])
            coeff = witness * B.inv()
            assert any(x for row in witness.tolist() for x in row)
            assert all(x.denominator == 1 for row in coeff.tolist() for x in row)
            norm = F(int((witness * G * witness.transpose())[0, 0]), p**2)
            assert norm == F(a["witness_squared_norm"])
            # Exact density upper bound from one actual vector, not an
            # unproved shortest-vector lower bound from lattice reduction.
            ratio2 = index**2 * (norm / 10) ** n
            assert ratio2 < 1 and a["rejected"]
            bounds.append(norm)
        summaries.append(
            {
                "dimension": n,
                "overlattice_index": index,
                "rejected_candidates": len(trials),
                "largest_supplied_minimum_upper_bound": str(max(bounds)),
                "scope": "Only these 24 candidates; no exclusion of the full family.",
            }
        )
    with np.load(P / "adaptive54.npz", allow_pickle=False) as data:
        rs, ws, ns = [
            data[k]
            for k in ["added_numerators", "witness_numerators", "norm_numerators"]
        ]
    assert rs.shape == ws.shape == (300, 54) and ns.shape == (300,)
    q3 = np.array([int(x) for x in (3 * q).tolist()[0]], dtype=np.int64)
    assert all(np.issubdtype(a.dtype, np.integer) for a in (rs, ws, ns))
    assert np.all((0 <= rs) & (rs < 3))
    assert np.max(np.abs(ws)) < 10**6 and np.max(np.abs(q3)) < 10**6
    for r, w, saved_norm in zip(rs, ws, ns):
        assert nmod_mat([r.tolist(), q3.tolist()], 3).rank() == 2
        rf, wf = fmpz_mat([r.tolist()]), fmpz_mat([w.tolist()])
        assert all(int(x) % 3 == 0 for row in (rf * G54).tolist() for x in row)
        # Membership in L + ZQ + Z(r/3) is exactly this residue condition.
        assert any(
            np.all((w - a * r - b * q3) % 3 == 0) for a in range(3) for b in range(3)
        )
        norm = F(int((wf * G54 * wf.transpose())[0, 0]), 9)
        assert norm == F(int(saved_norm), 3) > 0
        assert 81 * (norm / 10) ** 54 < 1
    assert int(ns.max()) == 22
    return {
        "failed_55D_deformation_squared_norm": bad["exact_squared_norm"],
        "overlattice_screens": summaries,
        "adaptive54_screen": {
            "seed": 545600,
            "rejected_candidates": 300,
            "largest_supplied_minimum_upper_bound": "22/3",
            "BKZ_block_size": 32,
            "scope": "Only these 300 candidates; no exclusion of the full family.",
        },
        "all_counterexamples_checked_exactly": True,
    }


def audit():
    source, (G, q, _, _) = frontier.parents54_55(return_data=True)
    assert source["counterexample55"]["squared_norm"] == 6
    assert G.det() == 3**25
    D, alpha = F(28, 3), F(1, 21)
    s2 = D * alpha**2
    h2 = F(10, 9) - s2
    assert s2 == F(4, 189) and h2 == F(206, 189)
    # Q+L and 2Q+L both have minimum D. Interpolation between w and
    # w+Q therefore applies although Q itself is not in L.
    layer1 = D - D * alpha * (1 - alpha) + h2
    # Reverse triangle inequality, using ||2 alpha Q||=2 alpha sqrt(D).
    layer2 = D * (1 - 2 * alpha) ** 2 + 4 * h2
    layer3 = 9 * (s2 + h2)
    assert layer1 == layer3 == 10 and layer2 > 10
    assert 0 < alpha < F(1, 2) and 9 * s2 <= F(10, 4)
    assert 16 * h2 > 10
    shift = q * 22 / 21
    cross = shift * fmpq_mat(G)
    bottom = (shift * fmpq_mat(G) * shift.transpose())[0, 0] + fmpq_mat([[str(h2)]])[
        0, 0
    ]
    gram = build(G, cross.tolist()[0], bottom)
    assert F(str(gram.det())) == 3**25 * h2
    witness = fmpq_mat([[-int(x) for x in (3 * q).tolist()[0]] + [3]])
    assert witness * gram * witness.transpose() == fmpq_mat([[10]])
    density = F(5, 2) ** 55 / (3**25 * h2)
    previous = json.loads((P.parent / "tilted55_20261116/results.json").read_text())
    prev2 = F(previous["density_squared"])
    repaired = F(previous["corrected_untilted_lamination_density_squared"])
    printed = F(previous["reported_external_claim_density_squared"])
    assert density / prev2 == F(1547, 1545)
    assert density / repaired == F(105, 103)
    assert density / printed == F(63, 103) < 1
    return {
        "status": "verified_cross_coset_tilt_55D",
        "dimension": 55,
        "name": "Cross-coset tilted lamination of Leterrier54",
        "kind": "new_explicit_lattice_source_comparison_disputed",
        "minimum_squared": 10,
        "center_density": decsqrt(density),
        "density_squared": str(density),
        "gram_determinant": str(gram.det()),
        "shift": "(22/21) Q",
        "vertical_height_squared": str(h2),
        "layer_minimum_squared_bounds": {
            "1": str(layer1),
            "2": str(layer2),
            "3": str(layer3),
            "all_absolute_layers_at_least_4": str(16 * h2),
        },
        "minimum_attainment_vector": list(map(int, witness.tolist()[0])),
        "previous_55D_density_squared": str(prev2),
        "ratio_squared_to_previous": "1547/1545",
        "gain_over_previous": gain(density / prev2),
        "corrected_untilted_lamination_density_squared": str(repaired),
        "ratio_squared_to_repaired_lamination": "105/103",
        "gain_over_repaired_lamination": gain(density / repaired),
        "reported_external_claim_density_squared": str(printed),
        "ratio_squared_to_printed_claim": "63/103",
        "world_record_claimed": False,
        "inherited_inputs": "Leterrier54 minimum10 and the minima28/3 of both nontrivial non-narrow cosets; the source-to-Gram reconstruction is repeated.",
        "maximum_parent_or_target_dimension": 56,
        "all_layer_and_density_checks_exact": True,
        "negative_results": negative_audits(G, q),
    }, gram


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-result", action="store_true")
    args = ap.parse_args()
    result, gram = audit()
    matrix = [[str(x) for x in row] for row in gram.tolist()]
    if args.write_result:
        (P / "results.json").write_text(json.dumps(result, indent=2) + "\n")
        (P / "gram55.json").write_text(json.dumps(matrix, indent=2) + "\n")
    else:
        assert result == json.loads((P / "results.json").read_text())
        assert matrix == json.loads((P / "gram55.json").read_text())
    print("55D center density", result["center_density"])
    print("Gain over previous:", result["gain_over_previous"]["percent"], "%")
    print("348 overlattices and one 55D deformation rejected exactly; no new SOTA.")
