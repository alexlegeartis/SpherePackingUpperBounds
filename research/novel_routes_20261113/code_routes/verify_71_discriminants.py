#!/usr/bin/env python3
"""Exact small-discriminant Gauss sums for the71D determinant corollary."""
from pathlib import Path
from math import gcd
import json
import sympy as sp


def phase(r):
    return sp.expand_complex(sp.exp(sp.pi * sp.I * r))


def brown(sum_value, cardinality):
    matches = [s for s in range(8)
               if sp.simplify(sum_value - sp.sqrt(cardinality) * phase(sp.Rational(s, 4))) == 0]
    assert len(matches) == 1
    return matches[0]


def main():
    cyclic = []
    for d in (2, 4, 6):
        for a in range(1, 2*d):
            if gcd(a, d) != 1:
                continue
            value = sum(phase(sp.Rational(a*x*x, d)) for x in range(d))
            cyclic.append({"order": d, "generator_norm_numerator": a,
                           "brown_invariant": brown(value, d)})
        matching = [x["generator_norm_numerator"] for x in cyclic
                    if x["order"] == d and x["brown_invariant"] == 7]
        assert matching == [2*d - 1]

    elementary = []
    for a in range(4):
        for b in range(4):
            for c in range(2):
                if (a*b - c*c) % 2 == 0:
                    continue
                value = sum(phase(sp.Rational(a*x+b*y, 2) + c*x*y)
                            for x in range(2) for y in range(2))
                s = brown(value, 4)
                assert s % 2 == 0
                elementary.append({"a": a, "b": b, "c": c, "brown_invariant": s})

    result = {
        "cyclic_forms": cyclic,
        "elementary_order4_forms": elementary,
        "unique_order_d_brown7_form_is_opposite_rank1_d_for_d2_4_6": True,
        "no_elementary_order4_brown7_form": True,
    }
    path = Path(__file__).with_name("discriminants71_verified.json")
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
