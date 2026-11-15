"""Exact 84D/85D comparisons to Chen (2013), Theorem 2.4, Section 3.

This verifies gains over a named published guarantee. It does not certify
exhaustiveness of the literature search or the unknown actual minimum of
Chen's lattice. Only the standard library is required; never run with -O.
"""
from decimal import Decimal, localcontext
from fractions import Fraction as F
from pathlib import Path
import argparse
import json

if not __debug__:
    raise RuntimeError("Run with assertions enabled.")

P = Path(__file__).resolve().parent


def sqrt_decimal(x):
    with localcontext() as ctx:
        ctx.prec = 50
        return (Decimal(x.numerator)/Decimal(x.denominator)).sqrt()


def audit():
    ours = json.loads((P/"counting80_results.json").read_text())
    assert ours["all_arithmetic_exact"] and ours["existential_not_explicit"]
    rows = []
    for n in (84, 85):
        m, ell, k = 10, 89, 1
        assert m < F(n, 2) and ell >= n+1
        assert all(ell % d for d in range(2, 10))  # primality of 89
        # The code used by Chen has just two words: zero and 84 ones
        # (and, for n=85, one zero). Its distance is therefore exactly 84.
        word = [1]*84+[0]*(n-84)
        assert len(word) == n and sum(word) == 84 >= 8*m
        assert sum(word) % 2 == 0
        # Independent evaluations of Theorem 2.4 and of radius/covolume.
        theorem2 = F(2)**(2*k-n)*F(m)**n/F(ell**(2*m-2)*(n+1))
        squared_radius = F(8*m, 4)
        squared_covolume = 2**(2*(n-k))*ell**(2*m-2)*(n+1)
        baseline2 = squared_radius**n/squared_covolume
        assert baseline2 == theorem2 == F(4*5**n, 89**18*(n+1))
        row = next(r for r in ours["rows"] if r["dimension"] == n)
        our2 = F(row["center_density_squared_exact"])
        assert our2 == (F(2**90, 625) if n == 84 else F(2**94, 9375))
        assert our2 == F(2**n, 1)/F(row["gram_determinant_exact"])
        ratio2 = our2/baseline2
        lo, hi = {84: ("899.675712", "899.675714"),
                  85: ("364.438552", "364.438554")}[n]
        assert (1+F(lo)/100)**2 < ratio2 < (1+F(hi)/100)**2
        with localcontext() as ctx:
            ctx.prec = 45
            gain = 100*(sqrt_decimal(ratio2)-1)
            logold = sqrt_decimal(baseline2).ln()/Decimal(2).ln()
            lognew = sqrt_decimal(our2).ln()/Decimal(2).ln()
        rows.append({
            "dimension": n,
            "comparison_source": "Hao Chen (2013), DOI 10.5802/jtnb.825, Theorem 2.4 and Section 3",
            "comparison_kind": "published_guarantee_reconstructed_from_theorem_and_specified_code",
            "comparison_center_density_squared_exact": str(baseline2),
            "comparison_center_density": str(sqrt_decimal(baseline2)),
            "comparison_log2_center_density": str(logold),
            "paper_table3_1_log2_display": "37.006" if n == 84 else "38.1616",
            "paper_display_note": "The displayed logarithm differs from evaluating the stated theorem; the exact theorem formula is used here.",
            "code_parameters": [n, k, 84],
            "generalized_craig_parameters": {"m": m, "ell": ell},
            "guaranteed_squared_minimum": 8*m,
            "actual_old_lattice_minimum_computed_here": False,
            "our_center_density_squared_exact": str(our2),
            "our_center_density": str(sqrt_decimal(our2)),
            "our_log2_center_density": str(lognew),
            "ratio_squared_exact": str(ratio2),
            "ratio": str(sqrt_decimal(ratio2)),
            "increase_percent": str(gain),
            "increase_percent_interval_verified_exactly": [lo, hi],
            "our_packing_explicit": False,
            "status": "verified_gain_over_named_published_guarantee",
            "literature_search_found_later_stronger_same_dimension_guarantee": False,
            "world_record_confirmed": False,
        })
    return {
        "checked_date": "2026-11-15",
        "source_doi": "https://doi.org/10.5802/jtnb.825",
        "source_pdf": "https://www.numdam.org/item/JTNB_2013__25_1_59_0.pdf",
        "source_locator": "Journal de Theorie des Nombres de Bordeaux 25 (2013), 59-70; Theorem 2.4 p.64, Section 3 p.65, Table 3.1 p.66",
        "all_gain_arithmetic_exact": True,
        "literature_review": "LITERATURE_STATUS.md",
        "source_manifest": "literature_sources.json",
        "rows": rows,
        "global_literature_search_exhaustive": False,
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-result", action="store_true")
    args = ap.parse_args()
    result = audit()
    path = P/"external_comparisons.json"
    if args.write_result:
        path.write_text(json.dumps(result, indent=2)+"\n")
    else:
        assert json.loads(path.read_text()) == result, "Comparison report is stale."
    for row in result["rows"]:
        print(f'{row["dimension"]}D: Chen guarantee={row["comparison_center_density"]}; '
              f'ours/Chen={row["ratio"]}; gain={row["increase_percent"]}%; '
              f'world_record_confirmed={row["world_record_confirmed"]}')
