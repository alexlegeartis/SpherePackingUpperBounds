#!/usr/bin/env python3
"""Exact rational torus-hole bounds propagated from Schuett's56D lattice."""

from fractions import Fraction as F
from math import factorial, gamma, pi, sqrt
from pathlib import Path
import json


HERE = Path(__file__).resolve().parent


def atan_partial(q, terms):
    return sum((F((-1)**k, (2*k + 1) * q**(2*k + 1))
                for k in range(terms)), F(0))


def main():
    # Machin's identity and alternating-series bounds certify the pi bound.
    pi_upper = 16*atan_partial(5, 15) - 4*atan_partial(239, 14)
    assert pi_upper < F(355, 113)
    parent = F(5**28, 2**28 * 29**4)
    exact = []
    for n, s, R2, same_color_ratio in [
        (57, 2, F(2327, 1000), 4),
        (58, 3, F(227, 100), 3),
    ]:
        mass_upper = (s - 1) * parent * F(355, 113)**28 * R2**28 / factorial(28)
        assert mass_upper < 1
        t2 = 4 - R2
        assert same_color_ratio * t2 >= 4
        squared_density = parent**2 / t2 if n == 57 else 4 * parent**2 / (3 * t2**2)
        exact.append({
            "dimension": n,
            "colors": s,
            "torus_separation_squared": str(R2),
            "vertical_minimum_squared": str(t2),
            "union_volume_upper_bound_exact": str(mass_upper),
            "union_volume_upper_bound_approx": float(mass_upper),
            "same_color_squared_minimum_ratio": same_color_ratio,
            "center_density_squared_exact": str(squared_density),
            "center_density_approx": sqrt(float(squared_density)),
            "certified_from_published_parent": True,
            "historical_priority_established": False,
        })

    volume = pi**28 / factorial(28)
    optimum_R2_s2 = (1 / (float(parent) * volume))**F(1, 28)
    optimum_R2_s3 = (1 / (2 * float(parent) * volume))**F(1, 28)
    limiting = [
        {"dimension": 57, "R2": optimum_R2_s2,
         "center_density_approx": float(parent) / sqrt(4-optimum_R2_s2)},
        {"dimension": 58, "R2": optimum_R2_s3,
         "center_density_approx": 2*float(parent)/(sqrt(3)*(4-optimum_R2_s3))},
    ]
    result = {
        "parent_dimension": 56,
        "parent_center_density_exact": str(parent),
        "parent_center_density_approx": float(parent),
        "primary_table": "https://cohn.mit.edu/sphere-packing/",
        "original_construction": "Schuett-Shioda, Mordell-Weil Lattices (2019), p400",
        "exact_algebraic_bounds": exact,
        "limiting_formula_approximations": limiting,
        "scope": "Derived comparison bounds; no claim these dimension-specific consequences were separately published.",
    }
    (HERE / "torus_layer_bounds_verified.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
