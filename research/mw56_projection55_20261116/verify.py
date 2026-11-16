"""Exact witness/coloring audit for 24 projection families from published MW56.

The full printed Leterrier55 number is the target, not a proved parent input.
The upper bounds here apply only to the specified translate families.
"""

from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext
from math import gcd
from functools import reduce
import argparse
import importlib.util
import json

import numpy as np
from flint import fmpz_mat

if not __debug__:
    raise RuntimeError("Run without -O.")
P = Path(__file__).resolve().parent


def decsqrt(x):
    with localcontext() as ctx:
        ctx.prec = 65
        return str((Decimal(x.numerator) / x.denominator).sqrt())


def bits(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def audit():
    spec = importlib.util.spec_from_file_location(
        "mw56_primary_source", P.parent / "cm58_20261116/verify.py"
    )
    cm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cm)
    G, _, _, _ = cm.parent()
    assert G.det() == 29**8
    H = 29 * G.inv()
    assert all(x.denominator == 1 for row in H.tolist() for x in row)
    H = fmpz_mat([[int(x) for x in row] for row in H.tolist()])
    target_data = json.loads((P / "targets.json").read_text())
    targets = {a["dimension"]: F(a["density_squared"]) for a in target_data["rows"]}
    assert targets == {
        54: F(5, 2) ** 54 / 3**25,
        55: F(5, 2) ** 55 / (2 * 3**24),
        56: F(5, 2) ** 56 / 29**8,
    }
    assert not target_data["rows"][1]["verified_lower_bound"]
    old = json.loads((P.parent / "tilted55_20261116/results.json").read_text())
    assert F(old["reported_external_claim_density_squared"]) == targets[55]
    target = targets[55]
    with np.load(P / "certificate.npz", allow_pickle=False) as data:
        saved = {k: data[k] for k in data.files}
    assert all(np.issubdtype(a.dtype, np.integer) for a in saved.values())
    normals, numerators = saved["normals"], saved["normal_numerators"]
    assert normals.shape == (24, 56) and numerators.shape == (24,)
    colorings = json.loads((P / "colorings.json").read_text())
    assert len(colorings) == 24
    rows, total_witnesses = [], 0
    for case, (normal, numerator, cert) in enumerate(
        zip(normals, numerators, colorings)
    ):
        assert case == cert["case"]
        b = fmpz_mat([list(map(int, normal))])
        assert reduce(gcd, map(abs, b.tolist()[0])) == 1
        m = int(numerator)
        assert int((b * H * b.transpose())[0, 0]) == m
        v0 = b * H
        assert reduce(gcd, map(abs, v0.tolist()[0])) == 1
        assert (v0 * G * v0.transpose())[0, 0] == 29 * m
        assert (v0 * b.transpose())[0, 0] == m
        col, T = b.transpose().hnf(transform=True)
        assert abs(T.det()) == 1 and col == T * b.transpose()
        assert col.tolist() == [[1]] + [[0]] * 55
        N = fmpz_mat(T.tolist()[1:])
        Q = N * G * N.transpose()
        assert N * b.transpose() == fmpz_mat(55, 1)
        assert Q.det() == 29**7 * m
        W = saved[f"witnesses_{case:02d}"]
        upper = list(map(int, saved[f"distance_numerators_{case:02d}"]))
        assert W.shape == (m, 56) and len(upper) == m
        assert upper[0] == 10 * m
        assert all(upper[j] == upper[-j] > 0 for j in range(1, m))
        for j, w in enumerate(W):
            wf = fmpz_mat([list(map(int, w))])
            label = int((wf * b.transpose())[0, 0])
            norm_numerator = m * int((wf * G * wf.transpose())[0, 0]) - 29 * label**2
            assert label % m == j and norm_numerator == upper[j] > 0
        total_witnesses += m
        # Exact upper bounds at all critical distance thresholds. These
        # cover arbitrary positive packing minima, not just minimum10.
        critical = sorted({u for u in upper if 0 < u <= 10 * m})
        assert [a["cutoff_numerator"] for a in cert["thresholds"]] == critical
        maximum = F(0)
        maximizing = None
        for item in cert["thresholds"]:
            cutoff = item["cutoff_numerator"]
            allowed = [j for j in range(1, m) if upper[j] >= cutoff]
            mask = sum(1 << j for j in allowed)
            adjacency = [sum(1 << ((v + j) % m) for j in allowed) for v in range(m)]
            colors = list(map(int, item["color_classes"]))
            used = 0
            for color in colors:
                assert color > 0 and color < 1 << m and color & ~mask == 0
                assert used & color == 0
                used |= color
                assert all(adjacency[v] & color == 0 for v in bits(color))
            assert used == mask and item["clique_upper"] == 1 + len(colors)
            count = 1 + len(colors)
            density2 = count**2 * F(cutoff, 4 * m) ** 55 / (29**7 * m)
            if density2 > maximum:
                maximum = density2
                maximizing = {
                    "minimum_squared_upper": str(F(cutoff, m)),
                    "translate_count_upper": count,
                }
        assert maximum < target
        rows.append(
            {
                "case": case,
                "normal_norm_squared": str(F(m, 29)),
                "quotient_order": m,
                "period_gram_determinant": str(29**7 * m),
                "all_translate_subsets_excluded_above_target": True,
                "density_upper_squared": str(maximum),
                "center_density_upper": decsqrt(maximum),
                "maximizing_upper_bound": maximizing,
                "thresholds_checked": len(critical),
            }
        )
    # These genuinely exceeded the printed target on the minimum10
    # hypothesis. Supplied short vectors reject that hypothesis exactly.
    proposals = []
    for case, count, label in [
        (1, 2, 28),
        (5, 2, 52),
        (8, 3, 36),
        (14, 2, 4),
        (17, 2, 60),
    ]:
        m = int(numerators[case])
        nominal = count**2 * F(5, 2) ** 55 / (29**7 * m)
        obstruction = F(int(saved[f"distance_numerators_{case:02d}"][label]), m)
        assert nominal > target and obstruction < 10
        assert nominal * (obstruction / 10) ** 55 < target
        proposals.append(
            {
                "case": case,
                "translates": count,
                "density_if_minimum10": decsqrt(nominal),
                "nominal_ratio_squared_to_printed_target": str(nominal / target),
                "verified_short_vector_squared_norm": str(obstruction),
                "status": "rejected; the displayed nominal density is not a lower bound",
            }
        )
    ceiling = max(F(a["density_upper_squared"]) for a in rows)
    return {
        "status": "verified_finite_family_exclusion_against_full_printed_target",
        "parent": "Published MW56, minimum10 determinant29^8; source-to-Gram audit repeated",
        "dimension": 55,
        "printed_target_squared": str(target),
        "printed_target_center_density": decsqrt(target),
        "printed_target_used_as_proved_parent": False,
        "projection_families_checked": len(rows),
        "individual_distance_witnesses_checked": total_witnesses,
        "all_family_upper_bounds_below_target": True,
        "largest_family_density_upper_squared": str(ceiling),
        "largest_family_center_density_upper": decsqrt(ceiling),
        "rows": rows,
        "rejected_record_sized_proposals": proposals,
        "scope": "Only the 24 supplied normal directions and subsets of their projected MW56 cosets, with the section lattice as period. Other directions, shifts and constructions remain open.",
        "new_lower_bound_claimed": False,
        "maximum_parent_or_target_dimension": 56,
        "all_proof_checks_exact": True,
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-result", action="store_true")
    args = ap.parse_args()
    result = audit()
    if args.write_result:
        (P / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    else:
        assert result == json.loads((P / "results.json").read_text())
    print("Printed 55D target:", result["printed_target_center_density"])
    print(
        "24 MW56 projection families excluded; exact distance witnesses:",
        result["individual_distance_witnesses_checked"],
    )
    print(
        "Largest family density upper:", result["largest_family_center_density_upper"]
    )
