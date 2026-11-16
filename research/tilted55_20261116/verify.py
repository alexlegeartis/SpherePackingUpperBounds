"""Exact 55D tilted lamination and the short-period correction.

No world-record claim: the published 116783 claim is retained as a disputed
claim, not silently replaced by a weaker comparison.
"""

from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext
import argparse
import json
import sys

import numpy as np
from flint import fmpz_mat, fmpq_mat
from flint import nmod_mat

if not __debug__:
    raise RuntimeError("Run without -O.")
P = Path(__file__).resolve().parent
sys.path.insert(0, str(P.parent / "cm58_20261116"))
import verify_frontier as frontier


def sqrt_decimal(x):
    with localcontext() as ctx:
        ctx.prec = 70
        return str((Decimal(x.numerator) / Decimal(x.denominator)).sqrt())


def audit():
    report, (G, q, roots, _) = frontier.parents54_55(return_data=True)
    assert report["counterexample55"]["squared_norm"] == 6
    screen = json.loads((P / "parity_screen.json").read_text())
    rejected = next(a for a in screen["rows"] if a["dimension"] == 56)["candidates"]
    assert len(rejected) == 3
    with np.load(
        P.parent / "cm58_20261116/certificate.npz", allow_pickle=False
    ) as data:
        G56 = fmpz_mat(data["G"].tolist())
        seed56 = data["point_coordinates"][0]
    assert G56.det() == 29**8
    maps56 = np.load(
        P.parent / "cm58_20261116/frontier/isometries56.npz", allow_pickle=False
    )["maps"]
    known56, _ = frontier.orbit(G56, seed56, maps56, 10)
    assert len(known56) == 2436
    for candidate in rejected:
        checks = candidate["checks"]
        witness = fmpz_mat([candidate["counterexample"]])
        assert nmod_mat(checks, 2).rank() == 7
        assert np.all(np.any((known56 @ np.array(checks).T) % 2, axis=1))
        assert all(
            int(x) % 2 == 0
            for row in (witness * fmpz_mat(checks).transpose()).tolist()
            for x in row
        )
        assert (witness * G56 * witness.transpose())[0, 0] == 10
        assert candidate["found_norm_upper"] == 10 and candidate["rejected"]
    assert G.det() == 3**25 and q.nrows() == 1 and q.ncols() == 54
    q3 = 3 * q
    assert all(x.denominator == 1 for row in q3.tolist() for x in row)
    v = fmpq_mat([roots.tolist()[0]])
    assert (v * fmpq_mat(G) * v.transpose())[0, 0] == 10
    assert F(str((q * fmpq_mat(G) * q.transpose())[0, 0])) == F(28, 3)
    alpha = F(2, 45)
    s = alpha.numerator * v / alpha.denominator
    s2 = 10 * alpha**2
    h2 = F(10, 9) - s2
    assert h2 == F(442, 405) > 0
    D = F(28, 3)
    # For 0<=t<=1 and w,w+v in the same coset, the exact identity
    # ||w+t v||²=(1-t)||w||²+t||w+v||²-t(1-t)||v||²
    # gives D-10t(1-t), using only the inherited coset lower bound.
    bounds = {
        1: D - 10 * alpha * (1 - alpha) + h2,
        2: D - 10 * (2 * alpha) * (1 - 2 * alpha) + 4 * h2,
        3: 9 * (s2 + h2),
        4: 16 * h2,
    }
    assert bounds[1] == bounds[3] == 10 and bounds[2] > 10 and bounds[4] > 10
    assert 0 < 2 * alpha < 1 and 2 * (3 * alpha) <= 1
    # For layer3, ||3s|| <= sqrt(10)/2 guarantees that zero is a nearest
    # parent point. Every layer |k|>=4 has vertical squared norm>=16h².
    shift = q + s
    label = fmpq_mat(G) * shift.transpose()
    full = fmpq_mat(55, 55)
    for i in range(54):
        for j in range(54):
            full[i, j] = G[i, j]
        full[i, 54] = full[54, i] = label[i, 0]
    full[54, 54] = (shift * fmpq_mat(G) * shift.transpose())[0, 0] + fmpq_mat(
        [[str(h2)]]
    )[0, 0]
    assert F(str(full.det())) == 3**25 * h2
    w = fmpq_mat([[-int(x) for x in q3.tolist()[0]] + [3]])
    assert (w * full * w.transpose())[0, 0] == 10
    d54 = F(5, 2) ** 54 / 3**25
    corrected = d54 * F(9, 4)
    reported = d54 * F(15, 4)
    density = F(5, 2) ** 55 / F(str(full.det()))
    assert density / corrected == F(225, 221) > 1
    assert density < reported
    # The original unsheared construction has the vertical norm6 vector.
    # For its genuine minimum6, all other layers satisfy the required
    # minimum using D and the inherited parent minimum10.
    assert D + F(2, 3) == 10 and 9 * F(2, 3) == 6
    old_actual = reported * F(3, 5) ** 55
    with localcontext() as ctx:
        ctx.prec = 70
        gain = 100 * (Decimal(sqrt_decimal(F(225, 221))) - 1)
        lo = F(int(gain * 10**6), 10**6)
        hi = lo + F(1, 10**6)
    assert (1 + lo / 100) ** 2 <= F(225, 221) < (1 + hi / 100) ** 2
    return {
        "status": "verified_new_55D_lattice_and_rejected_old_period_claim",
        "world_record_claimed": False,
        "record_scope": "Improves the repaired order-three lamination; does not exceed the printed 116783 claim. The displayed construction for that claim has a norm6 counterexample.",
        "dimension": 55,
        "name": "Tilted order-three lamination of Leterrier54",
        "kind": "new_explicit_lattice_source_comparison_disputed",
        "minimum_squared": 10,
        "density_squared": str(density),
        "center_density": sqrt_decimal(density),
        "gram_determinant": str(full.det()),
        "alpha": str(alpha),
        "vertical_height_squared": str(h2),
        "layer_minimum_squared_bounds": {str(k): str(x) for k, x in bounds.items()},
        "minimum_attainment_vector": list(map(int, w.tolist()[0])),
        "corrected_untilted_lamination_density_squared": str(corrected),
        "corrected_untilted_lamination_center_density": sqrt_decimal(corrected),
        "ratio_squared_to_corrected_lamination": "225/221",
        "gain_percent_interval": [str(lo), str(hi)],
        "reported_external_claim_density_squared": str(reported),
        "reported_external_claim_center_density": sqrt_decimal(reported),
        "reported_external_claim_used_as_valid_comparator": False,
        "rejected_previous_55D_gram_actual_minimum_squared": 6,
        "rejected_previous_55D_gram_actual_center_density": sqrt_decimal(old_actual),
        "counterexample_to_previous_minimum10_claim": report["counterexample55"],
        "inherited_inputs": "Published Leterrier54 minimum10 and non-narrow coset minimum28/3; explicit source-to-Gram reconstruction is repeated.",
        "maximum_parent_or_target_dimension": 55,
        "all_layer_and_density_checks_exact": True,
    }, full


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-result", action="store_true")
    args = ap.parse_args()
    result, gram = audit()
    matrix = [[str(x) for x in r] for r in gram.tolist()]
    if args.write_result:
        (P / "results.json").write_text(json.dumps(result, indent=2) + "\n")
        (P / "gram55.json").write_text(json.dumps(matrix, indent=2) + "\n")
    else:
        assert result == json.loads((P / "results.json").read_text())
        assert matrix == json.loads((P / "gram55.json").read_text())
    print("55D minimum10, center density", result["center_density"])
    print("Old minimum10 claim rejected by norm6 vector; world-record claim withheld.")
