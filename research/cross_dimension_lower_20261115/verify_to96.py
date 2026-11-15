"""Counting extensions and classical controls, with a hard dimension cap of 96.

Standard-library exact checks; mathematical inputs/proofs are in TO96.md.
Classical controls derived here are not attributed as reported records.
"""
from decimal import Decimal, localcontext
from fractions import Fraction as F
from math import factorial, ceil, exp, lgamma, pi
from pathlib import Path
import argparse
import json
from verify_counting import theta, determinant, root_gram, decsqrt

if not __debug__:
    raise RuntimeError("Assertions must be enabled.")

P = Path(__file__).resolve().parent
MAX_DIMENSION = 96


def exceptional(r):
    assert r in (6, 7, 8)
    g = [[2*int(i == j) for j in range(r)] for i in range(r)]
    for i, j in [(i, i+1) for i in range(r-2)]+[(2, r-1)]:
        g[i][j] = g[j][i] = -1
    assert determinant(g) == {6: 3, 7: 2, 8: 1}[r]
    return g


def direct_sum(a, b):
    n, m = len(a), len(b)
    return [row+[0]*m for row in a]+[[0]*n+row for row in b]


def atan_interval(a, terms):
    s = sum((F((-1)**j, (2*j+1)*a**(2*j+1)) for j in range(terms)), F(0))
    t = F((-1)**terms, (2*terms+1)*a**(2*terms+1))
    return min(s, s+t), max(s, s+t)


def pi_interval():
    a, b = atan_interval(5, 25), atan_interval(239, 8)
    return 16*a[0]-4*b[1], 16*a[1]-4*b[0]


def hermite_upper(d, pi_lower):
    # Floating arithmetic proposes a centesimal upper bound; only the exact
    # rational-power comparison below certifies it.
    c = F(ceil(100*(2/pi)*exp(2*lgamma(d/2+2)/d)), 100)
    if d % 2 == 0:
        bpow = (2/pi_lower)**d*factorial(d//2+1)**2
    else:
        j = d//2+2
        g = F(factorial(2*j), 4**j*factorial(j))
        bpow = F(2)**d*g*g/pi_lower**(d-1)
    assert c**d > bpow
    return c


def section_floor(d2, d, c):
    # delta_(d-1)^2 >= (delta_d^2)^((d-1)/d)/gamma_d.
    with localcontext() as ctx:
        ctx.prec = 70
        x = Decimal(d2.numerator)/Decimal(d2.denominator)
        candidate = int((x.ln()*Decimal(d-1)/Decimal(d)).exp()
                        /(Decimal(c.numerator)/Decimal(c.denominator)))
    threshold = d2**(d-1)
    while (c*candidate)**d > threshold:
        candidate -= 1
    while (c*(candidate+1))**d <= threshold:
        candidate += 1
    assert (c*candidate)**d <= threshold < (c*(candidate+1))**d
    return F(candidate)


def audit():
    t, _ = theta(80)
    bad = sum(t[n//2] for n in (8, 12, 16))//2
    population = 2**79-2**39-1
    counted = []
    for r in range(6, 11):
        g = exceptional(r) if r <= 8 else direct_sum(exceptional(8), root_gram(r-8))
        assert all(g[i][j] == g[j][i] for i in range(r) for j in range(r))
        assert all(g[i][i] == 2 for i in range(r))
        assert all(determinant([row[:j] for row in g[:j]]) > 0 for j in range(1, r+1))
        # Integral, even, positive definite Gram with a diagonal 2: min=2.
        expectation = F((2**r-1)*bad, population)
        assert expectation < 1
        det = F(3, 2)**r*determinant(g)
        d2 = F(2**(80+r), 1)/det
        counted.append({"dimension": 80+r, "rank": r,
                        "root_gram": g, "gram_determinant_exact": str(det),
                        "expected_bad_classes_upper_bound_exact": str(expectation),
                        "minimum_squared_exact_inheriting_parent": 8,
                        "center_density_squared_exact": str(d2),
                        "center_density": decsqrt(d2),
                        "explicit_horizontal_representatives": False,
                        "world_record_confirmed": False})
    assert F((2**11-1)*bad, population) > 1  # no rank-11 conclusion
    n86 = counted[0]
    ours = F(n86["center_density_squared_exact"])
    assert ours == F(2**92, 3**7)
    chen = F(4*5**86, 89**18*87)
    ratio2 = ours/chen
    with localcontext() as ctx:
        ctx.prec = 45
        ratio = (Decimal(ratio2.numerator)/Decimal(ratio2.denominator)).sqrt()
        gain = 100*(ratio-1)
    assert ratio2 > 1
    gain_interval = ["116.264339", "116.264340"]
    assert (1+F(gain_interval[0])/100)**2 < ratio2
    assert ratio2 < (1+F(gain_interval[1])/100)**2

    # Theorem 13.35 of the Schuett--Shioda 2017 draft, p=47:
    # rank=92, min=16, Gram determinant=47^14.
    mw92 = F(2**184, 47**14)
    assert mw92 == F(4**92, 47**14)
    lo, hi = pi_interval()
    assert F("3.141592653589793238462643383279") < lo < hi
    assert hi < F("3.141592653589793238462643383280")
    sections = []
    d2 = mw92
    for d in range(92, 84, -1):
        c = hermite_upper(d, lo)
        d2 = section_floor(d2, d, c)
        sections.append({"dimension": d-1,
                         "method": "Derived MW92 section using Blichfeldt's Hermite bound",
                         "hermite_upper_used": str(c),
                         "center_density_squared_lower_bound_exact": str(d2),
                         "center_density_lower_bound": decsqrt(d2),
                         "is_dimension_specific_result_reported_by_source": False,
                         "world_record_confirmed": False})
    products = []
    for r in range(1, 5):
        # Common minimum 16. Multiplication of center densities gives the
        # same scale-independent result, with root density^2=1/(2^r det R).
        d2 = mw92/F(2**r*determinant(root_gram(r)))
        products.append({"dimension": 92+r,
                         "method": "Derived MW92 orthogonal product with A1/A2/A3/D4",
                         "center_density_squared_exact": str(d2),
                         "center_density": decsqrt(d2),
                         "is_dimension_specific_result_reported_by_source": False,
                         "world_record_confirmed": False})
    # The 87--90D all-singular extensions are below unconditional classical
    # MW92 section controls, even without using recent rank-88 claims.
    for row in counted[1:]:
        section = next(s for s in sections if s["dimension"] == row["dimension"])
        assert F(row["center_density_squared_exact"]) < F(section["center_density_squared_lower_bound_exact"])
    assert max([r["dimension"] for r in counted+sections+products]) == MAX_DIMENSION
    return {"maximum_dimension": MAX_DIMENSION,
            "counting_extensions": counted,
            "86D_comparison": {"source": "Chen (2013), Theorem 2.4 and Section 3",
                               "comparison_center_density_squared_exact": str(chen),
                               "comparison_center_density": decsqrt(chen),
                               "ratio_squared_exact": str(ratio2),
                               "ratio": str(ratio), "increase_percent": str(gain),
                               "gain_percent_interval_verified_exactly": gain_interval,
                               "world_record_confirmed": False},
            "MW92": {"rank": 92, "minimum_squared": 16, "gram_determinant": str(47**14),
                      "center_density_squared_exact": str(mw92), "center_density": decsqrt(mw92)},
            "classical_sections": sections, "classical_products": products,
            "all_checks_exact": True,
            "conditional_88D_parent_used_as_proof_input": False}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-result", action="store_true")
    args = ap.parse_args()
    result = audit()
    path = P/"to96_results.json"
    if args.write_result:
        path.write_text(json.dumps(result, indent=2)+"\n")
    else:
        assert json.loads(path.read_text()) == result, "Stored to96 report is stale."
    print("86D:", result["counting_extensions"][0]["center_density"])
    print("Gain over Chen's guarantee:", result["86D_comparison"]["increase_percent"], "%")
    for row in result["classical_sections"]:
        print(row["dimension"], "MW92 section >=", row["center_density_lower_bound"])
    for row in result["classical_products"]:
        print(row["dimension"], "MW92 product =", row["center_density"])
