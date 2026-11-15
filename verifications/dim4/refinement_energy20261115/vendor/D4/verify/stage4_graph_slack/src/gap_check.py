"""Exact comparisons for the broad-antipodal nine-point gate.

No floating point. A nonzero exit means a failed assertion.
"""

from fractions import Fraction

RHO = Fraction(41, 49)
BETA = Fraction(7, 8)
CUTOFF = Fraction(223, 735)
ONE_THIRD = Fraction(1, 3)


def projection_maximum():
    return (RHO - BETA * BETA) / (1 - BETA * BETA)


def main():
    q = projection_maximum()
    gap = ONE_THIRD - q
    saving = Fraction(1358 - 1208, 8 * 50_000_000)
    if q != CUTOFF:
        raise SystemExit(f"projection maximum {q} is not {CUTOFF}")
    if gap != Fraction(22, 735):
        raise SystemExit(f"gap {gap} is not 22/735")
    if not gap > 0:
        raise SystemExit("cutoff is not strictly below 1/3")
    if not BETA > RHO:
        raise SystemExit("zero-projection cosine is not above the directional cutoff")
    if 151 * 8 != 1208:
        raise SystemExit("ordered capacity is not 1208")
    if saving != Fraction(3, 8_000_000):
        raise SystemExit(f"conditional saving {saving} is not 3/8000000")
    print(f"Q(7/8,7/8)={q}")
    print(f"gap={gap}")
    print(f"7/8-41/49={BETA - RHO}")
    print("ordered=1208")
    print(f"conditional_saving={saving}")


if __name__ == "__main__":
    main()
