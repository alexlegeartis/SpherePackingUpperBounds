#!/usr/bin/env python3
"""Exact arithmetic behind the 72D short-normal section obstruction.

The mathematical input is the classical harmonic-theta modularity theorem
and C[E4,E6], not a numerical minimum search. This script checks the required
cusp-space basis, its Fourier coefficients, and the coefficient relation.
See MODULAR_OBSTRUCTION.md for the complete proof and its scope.
"""
from fractions import Fraction
import json
from pathlib import Path


def mul(a, b, n=3):
    return [sum(a[i] * b[k-i] for i in range(k+1)) for k in range(n+1)]


def power(a, k):
    out = [1, 0, 0, 0]
    for _ in range(k):
        out = mul(out, a)
    return out


def sigma(n, k):
    return sum(d**k for d in range(1, n+1) if n % d == 0)


def main():
    e4 = [1] + [240*sigma(n, 3) for n in range(1, 4)]
    e6 = [1] + [-504*sigma(n, 5) for n in range(1, 4)]
    e4_cubed = power(e4, 3)
    e6_squared = power(e6, 2)
    delta = [Fraction(a-b, 1728) for a, b in zip(e4_cubed, e6_squared)]
    assert delta == [0, 1, -24, 252]

    # S_38 = Delta M_26, and the two monomials of weight26 span M_26.
    monomials_26 = [(a, b) for a in range(7) for b in range(5)
                    if 4*a+6*b == 26]
    assert monomials_26 == [(2, 3), (5, 1)]
    f1 = mul(delta, mul(power(e4, 5), e6))
    f2 = mul(power(delta, 2), mul(power(e4, 2), e6))
    assert f1 == [0, 1, 672, -51084]
    assert f2 == [0, 0, 1, -72]
    # Their distinct initial orders show independence; dimension is2.
    assert all(f[3]+72*f[2]+2700*f[1] == 0 for f in (f1, f2))

    # For x^2 in {2,4,6}, w^2<=10, and nonzero integral inner product,
    # H_w(x)=(x,w)^2 - x^2*w^2/72 is strictly positive.
    margins = {str(s): str(Fraction(1)-Fraction(s*10, 72))
               for s in (2,4,6)}
    assert min(Fraction(v) for v in margins.values()) == Fraction(1,6)
    out = {
        "arithmetic_passed": True,
        "basis_F1_q0_to_q3": [int(x) for x in f1],
        "basis_F2_q0_to_q3": [int(x) for x in f2],
        "relation": "T6 + 72*T4 + 2700*T2 = 0",
        "positive_margins_at_w_norm10": margins,
        "theorem": "Every nonextremal even unimodular72D lattice and every "
                   "nonzero lattice vector w with w^2<=10 have a nonzero "
                   "vector of norm<=6 in w-perp.",
        "world_density_record": False,
        "inputs_not_proved_by_script": [
            "Harmonic theta of degree2 for an even unimodular72D lattice is in S_38(SL2Z).",
            "The level-one modular-form ring is C[E4,E6], and S_k=Delta M_(k-12)."
        ],
    }
    p = Path(__file__).with_name("modular_obstruction_verified.json")
    p.write_text(json.dumps(out, indent=2)+"\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
