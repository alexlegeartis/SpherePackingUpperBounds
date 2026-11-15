"""Verify the comparison with Vardy's REPORTED 59D guarantee, using rationals.

The source prints a decimal prefix for a lower bound. It gives no upper bound
on the actual Elkies lattice density. This checker must not certify such a gain.
The source transcription/interpretation is an input; the arithmetic is exact.
"""
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path
import argparse
import hashlib
import json

if not __debug__:
    raise RuntimeError("Run without -O; exact proof checks require assertions.")

P = Path(__file__).resolve().parent
SOURCES = P.parents[1] / "research/novel_routes_20261113/audits/benchmark_sources.json"


def decimal(value):
    return Decimal(value.numerator) / Decimal(value.denominator)


def ratio_endpoint_check(delta2, exponent, percent, direction):
    """Compare sqrt(delta2)/2**exponent with 1+percent/100 exactly."""
    ratio = 1 + percent / 100
    assert delta2 > 0 and ratio > 0 and exponent > 0
    # All quantities are positive, so taking this even power is reversible.
    power = exponent.denominator
    left = delta2 ** power
    right = Fraction(2) ** (2 * exponent.numerator) * ratio ** (2 * power)
    assert (left > right) if direction == ">" else (left < right)
    return {
        "exponent": str(exponent), "percent": str(decimal(percent)),
        "ratio_rational": str(ratio), "direction": direction,
        "delta_squared_power": power,
        "power_of_two": 2 * exponent.numerator,
        "ratio_power": 2 * power,
        "passed": True,
    }


def audit():
    result_path = P / "result.json"
    result = json.loads(result_path.read_text())
    source = json.loads(SOURCES.read_text())["V99"]
    assert result["dimension"] == 59 and result["all_checks_exact"]
    delta2 = Fraction(result["center_density_squared_exact"])
    assert delta2 == Fraction(3**33 * 401**59, 2**296 * 5**2 * 13**61)
    observation = source["verified_59D_observation"]
    assert observation["printed_log2_guarantee_prefix"] == "17.75.."
    assert observation["wording_is_lower_bound"]
    assert observation["section_results_reference"] == "10"
    assert observation["section_results_origin"] == "N. D. Elkies, personal communication, February 1996"

    # Let B be the numerical GUARANTEE whose logarithm has the printed prefix.
    # This is an interval for log2(B), not for log2(delta(actual E59)).
    lower_log, upper_log = Fraction("17.75"), Fraction("17.76")
    lower_percent, upper_percent = Fraction("23.675981"), Fraction("24.536216")
    checks = [
        ratio_endpoint_check(delta2, upper_log, lower_percent, ">"),
        ratio_endpoint_check(delta2, lower_log, upper_percent, "<"),
    ]
    # The paper also gives an exact but weaker comparator: delta(J_n)=2^(n-42).
    j59 = Fraction(2**(59-42))
    assert j59 == 131072
    j_checks = [
        ratio_endpoint_check(delta2, Fraction(17), Fraction("109.444114"), ">"),
        ratio_endpoint_check(delta2, Fraction(17), Fraction("109.444115"), "<"),
    ]
    with localcontext() as ctx:
        ctx.prec = 85
        center_density = decimal(delta2).sqrt()
        endpoint_gain_low = 100 * (center_density / (Decimal(2) ** decimal(upper_log)) - 1)
        endpoint_gain_high = 100 * (center_density / (Decimal(2) ** decimal(lower_log)) - 1)
        j_gain = 100 * (center_density / Decimal(j59.numerator) - 1)
        display = {
            "endpoint_gain_lower_approx": format(endpoint_gain_low, ".45f"),
            "endpoint_gain_upper_approx": format(endpoint_gain_high, ".45f"),
            "J59_gain_approx": format(j_gain, ".45f"),
        }
    return {
        "status": "verified_gain_interval_over_reported_guarantee",
        "dimension": 59,
        "source_id": "V99",
        "source_url": source["url"],
        "source_locator": source["locator"],
        "source_transcription_is_external_input": True,
        "comparison_target": "The lower-bound value reported for E59, not the actual E59 density",
        "used_as_main_comparator": False,
        "stronger_comparison": "literature_section_comparison.json (Gamma72 perpendicular to classical kappa13)",
        "reported_guarantee_log2_interval": {
            "lower": str(lower_log), "lower_inclusive": True,
            "upper": str(upper_log), "upper_inclusive": False,
        },
        "gain_percent_interval": {
            "lower": "23.675981", "lower_inclusive": False,
            "upper": "24.536216", "upper_inclusive": False,
        },
        "exact_integer_power_checks": checks,
        "display_only_decimals": display,
        "actual_E59_density_known_exactly": False,
        "improvement_over_actual_E59_established": False,
        "current_record_improvement_established": False,
        "reason": "Vardy p. 278 reports only a lower bound; reference 10 on p. 283 is personal correspondence. The upper endpoint above bounds the reported guarantee, not the actual lattice density.",
        "weaker_exact_external_control": {
            "name": "Vardy J59", "center_density_exact": str(j59),
            "locator": "Section 4, p. 278, delta(J_n)=2^(n-42)",
            "gain_percent_exact_expression": "100*(sqrt(" + str(delta2 / j59**2) + ")-1)",
            "gain_percent_interval": ["109.444114", "109.444115"],
            "exact_integer_power_checks": j_checks,
            "used_as_main_comparator": False,
            "reason": "Vardy explicitly reports denser Elkies sections in the same paragraph.",
        },
        "construction_result_sha256": hashlib.sha256(result_path.read_bytes()).hexdigest(),
        "all_arithmetic_checks_exact": True,
        "scope": "Checks comparison arithmetic against manually verified source facts. Replay verify.py and verify_independent.py for the packing itself.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-result", action="store_true")
    args = parser.parse_args()
    report = audit()
    output = P / "external_comparison.json"
    if args.write_result:
        output.write_text(json.dumps(report, indent=2) + "\n")
    else:
        assert json.loads(output.read_text()) == report, "External comparison report is stale."
    print(json.dumps({k: report[k] for k in (
        "status", "gain_percent_interval", "improvement_over_actual_E59_established",
        "current_record_improvement_established", "all_arithmetic_checks_exact",
    )}, indent=2))
