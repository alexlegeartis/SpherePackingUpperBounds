#!/usr/bin/env python3
"""Independent exact comparison for the 71D rank-one antipode formula.

The optimization proof is in RANK1_AUDIT.md. This verifies its decisive
rational inequalities, and independently checks the actual norm10 normal
against the archived72D Gram. It does not certify world-record priority.
"""
from fractions import Fraction
import json
from pathlib import Path
import numpy as np


def main():
    p = Path(__file__).resolve().parent
    source = p.parent / "higher_parent"
    z = np.load(source / "gamma72_simplex_rank1.npz")
    G = z["G"].astype(object)
    B = z["B"].astype(object)
    S = z["S"].astype(object)
    assert B.shape == (1,72)
    assert int((B@G@B.T)[0,0]) == 10
    from math import gcd
    assert gcd(*[int(x) for x in B[0]]) == 1
    assert sorted(int(x) for x in S[:,0]) == [0,1]

    d2 = Fraction(2,5)*Fraction(79,40)**71
    rank1_at8_squared = Fraction(1,2)*Fraction(63,32)**71
    s_ge3_continuous_ceiling_squared = Fraction(1,4)*Fraction(71,36)**71
    bare_squared = Fraction(2)**68
    ratio8 = Fraction(4,5)*Fraction(316,315)**71
    ratio3 = Fraction(8,5)*Fraction(711,710)**71
    assert d2/rank1_at8_squared == ratio8 > 1
    assert d2/s_ge3_continuous_ceiling_squared == ratio3 > 1
    assert d2 > bare_squared

    out = {
        "all_exact_comparisons_passed": True,
        "norm10_normal_checked": True,
        "normal_primitive_checked": True,
        "center_density_squared": str(d2),
        "center_density_approx": float(d2)**0.5,
        "ratio_to_norm8_two_coset_approx": float(ratio8)**0.5,
        "ratio_to_bare_norm8_section_approx": float(d2/bare_squared)**0.5,
        "rank1_antipode_minimum8_formula_optimum": True,
        "world_record_established": False,
        "scope": "Formula using only the parent minimum8, all primitive "
                 "rank-one normals in an even unimodular72D lattice, "
                 "all finite distinct dual-label sets. Does not optimize "
                 "exact glue-fiber minima or arbitrary71D packings."
    }
    (p/"rank1_antipode_verified.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))


if __name__ == "__main__":
    main()
