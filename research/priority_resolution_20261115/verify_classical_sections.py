"""BCS88 block sections: exact classical comparisons in 84--87D.

The BCS lattice theorem and E8 norm-doubler are published inputs. Exact
Reed--Solomon checks bind the specified construction and its block projection.
"""
from pathlib import Path
from itertools import combinations
from fractions import Fraction as F
import argparse
import json
import sys

if not __debug__:
    raise RuntimeError("Assertions must be enabled.")
P = Path(__file__).resolve().parent
sys.path.insert(0, str(P.parent/"cross_dimension_lower_20261115"))
from verify_counting import determinant, root_gram, decsqrt


def mul(a, b):
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & 16:
            a ^= 19  # x^4+x+1
    return out


def power(a, n):
    out = 1
    for _ in range(n):
        out = mul(out, a)
    return out


def nonsingular(a):
    a = [row[:] for row in a]
    n = len(a)
    for i in range(n):
        p = next((j for j in range(i, n) if a[j][i]), None)
        if p is None:
            return False
        a[i], a[p] = a[p], a[i]
        inverse = power(a[i][i], 14)
        a[i] = [mul(x, inverse) for x in a[i]]
        for j in range(i+1, n):
            c = a[j][i]
            a[j] = [x ^ mul(c, y) for x, y in zip(a[j], a[i])]
    return True


def audit():
    assert all(power(a, 15) == 1 for a in range(1, 16))
    codes = []
    for k in (10, 8, 4):
        g = [[power(a, j) for a in range(11)] for j in range(k)]
        assert g[0] == [1]*11
        count = 0
        for cols in combinations(range(11), k):
            assert nonsingular([[row[c] for c in cols] for row in g])
            count += 1
        codes.append({"length": 11, "dimension_over_GF16": k,
                      "minimum_distance": 12-k,
                      "all_maximal_minors_checked": count,
                      "each_coordinate_projection_surjective": True})
    assert sum(c["all_maximal_minors_checked"] for c in codes) == 506
    # m=8, n=11, a=3 in BCS Eq.(15). Scale lengths by 2 to minimum 8.
    density88 = F(1, 16)**11*2**(4*(10+8+4))
    assert density88 == 2**44
    # Before scaling, covolume = (sqrt(2)/2)^88 / delta = 2^-88.
    assert F(1, 2**44)/density88 == F(1, 2**88)
    assert F(1, 2**88)*2**88 == 1  # after length scaling by 2
    ours80 = json.loads((P.parent/"cross_dimension_lower_20261115/counting80_results.json").read_text())
    to96 = json.loads((P.parent/"cross_dimension_lower_20261115/to96_results.json").read_text())
    claimed = json.loads((P.parent/"cross_dimension_lower_20261115/claimed88_comparisons.json").read_text())
    rows = []
    for k in range(1, 5):
        n = 88-k
        K = [[2*x for x in row] for row in root_gram(k)]
        det = determinant(K)
        d2 = F(2**n, 1)/det
        if n in (84, 85):
            our2 = F(next(r for r in ours80["rows"] if r["dimension"] == n)["center_density_squared_exact"])
            assert our2 > d2
        else:
            our2 = F(claimed["rows"][0 if n == 86 else 1]["center_density_squared_conditional"])
            assert our2 < d2
            if n == 86:
                assert F(to96["counting_extensions"][0]["center_density_squared_exact"]) < d2
        rows.append({"dimension": n, "normal_dual_rank": k,
                     "normal_dual_gram": K, "normal_dual_determinant": str(det),
                     "classical_center_density_squared_exact": str(d2),
                     "classical_center_density": decsqrt(d2),
                     "our_compared_center_density_squared": str(our2),
                     "our_to_classical_ratio_squared_exact": str(our2/d2),
                     "our_to_classical_ratio": decsqrt(our2/d2),
                     "our_bound_beats_classical_section": our2 > d2,
                     "numerical_dimension_specific_value_reported_in_source": False})
    return {"source": "Bos-Conway-Sloane (1982), Theorem 2, Eq.(15), E8 examples",
            "source_url": "https://neilsloane.com/doc/Me103.pdf",
            "codes": codes, "BCS88_center_density_exact": "2^44",
            "BCS88_minimum_after_length_scaling": 8,
            "BCS88_covolume_after_length_scaling": 1,
            "normal_block_input": "Dual intersection with any 8D coordinate block is sqrt(2)*E8, from surjectivity of all three code projections",
            "rows": rows, "all_arithmetic_checks_exact": True,
            "86D_and_87D_record_candidates_rejected": True,
            "Sakamoto88_minimum_claim_used": False}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-result", action="store_true")
    args = ap.parse_args()
    result = audit()
    path = P/"classical_bcs88_sections.json"
    if args.write_result:
        path.write_text(json.dumps(result, indent=2)+"\n")
    else:
        assert json.loads(path.read_text()) == result, "Stored comparison is stale."
    for row in result["rows"]:
        print(row["dimension"], row["classical_center_density"], "our/classical", row["our_to_classical_ratio"])
