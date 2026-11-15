#!/usr/bin/env python3
"""Independent exact proof of the moderate graph capacity 1880.

Only input dependency: moderate_exterior_graph_certificate.json.
Only runtime dependency: the Python standard library.
No archived verifier modules are imported or executed.

Reconstruct U_k by its explicit binomial formula, then use ordinary exact
rational Euclidean Sturm sequences. Positive integer rescaling of each
sequence entry preserves Sturm sign variations. All divisions and endpoint
evaluations are rational. A completed-square bound over REAL cardinalities
already suffices for the even ordered-edge capacity.
"""
import argparse
from fractions import Fraction as Q
from functools import reduce
import hashlib
import json
from math import comb, gcd, lcm
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def primitive(p):
    """Multiply by a strictly positive rational to get primitive integers."""
    p = trim(list(map(Q, p)))
    den = lcm(*(c.denominator for c in p))
    integers = [int(c * den) for c in p]
    content = reduce(gcd, (abs(c) for c in integers), 0)
    if content:
        integers = [c // content for c in integers]
    return integers


def remainder(a, b):
    a, b = trim(list(map(Q, a))), trim(list(map(Q, b)))
    require(b != [0], "Polynomial division by zero")
    while a != [0] and len(a) >= len(b):
        offset = len(a) - len(b)
        factor = a[-1] / b[-1]
        for j, coefficient in enumerate(b):
            a[offset + j] -= factor * coefficient
        a = trim(a)
    return a


def sturm(p):
    p = primitive(p)
    derivative = primitive([i * p[i] for i in range(1, len(p))])
    sequence = [p, derivative]
    while sequence[-1] != [0]:
        nxt = primitive([-c for c in remainder(sequence[-2], sequence[-1])])
        if nxt == [0]:
            break
        sequence.append(nxt)
    return sequence


def evaluate(p, x):
    value = Q(0)
    for c in reversed(p):
        value = value * x + c
    return value


def sign(value):
    return (value > 0) - (value < 0)


def variations(signs):
    nonzero = [s for s in signs if s]
    return sum(x != y for x, y in zip(nonzero, nonzero[1:]))


def prove_strict_positive(p, left, right):
    sequence = sturm(p)
    lsigns = [sign(evaluate(row, left)) for row in sequence]
    rsigns = [sign(evaluate(row, right)) for row in sequence]
    roots = variations(lsigns) - variations(rsigns)
    require(lsigns[0] == rsigns[0] == 1, "Required positive endpoint fails")
    require(roots == 0, "Required polynomial has an interior root")
    return {
        "closed_interval": [str(left), str(right)],
        "left_endpoint_value": str(evaluate(p, left)),
        "right_endpoint_value": str(evaluate(p, right)),
        "sturm_degrees": [len(row) - 1 for row in sequence],
        "left_signs": lsigns,
        "right_signs": rsigns,
        "left_variations": variations(lsigns),
        "right_variations": variations(rsigns),
        "interior_root_count": roots,
        "conclusion": "STRICTLY_POSITIVE_ON_ENTIRE_CLOSED_INTERVAL",
    }


def verify(path):
    raw = path.read_bytes()
    data = json.loads(raw)
    require(data["schema"] == "S3_SINGLE_COLOR_EXTERIOR_GRAPH_V1", "Wrong schema")
    require(data["ambient_dimension"] == 4, "Wrong ambient dimension")
    require(data["minimum_pair_distance_squared"] == "1", "Wrong separation")
    harmonic = list(map(Q, data["harmonic_coefficients"]))
    degree = data["harmonic_degree"]
    require(degree == 16 and len(harmonic) == 17, "Incomplete degree-16 coefficients")
    require(harmonic[0] > 0 and all(c >= 0 for c in harmonic), "Harmonic positivity fails")

    # U_k(t) = sum_{j=0}^{floor(k/2)} (-1)^j binom(k-j,j) (2t)^(k-2j).
    power = [Q(0)] * 17
    for k, coefficient in enumerate(harmonic):
        for j in range(k // 2 + 1):
            exponent = k - 2 * j
            power[exponent] += coefficient * Q((-1)**j * comb(k-j, j) * 2**exponent, k+1)
    require(power == list(map(Q, data["power_coefficients"])), "Harmonic/power identity fails")
    diagonal = evaluate(power, Q(1))
    require(diagonal == sum(harmonic), "Z_k(1)=1 identity fails")
    require(harmonic[0] == Q(109916026941, 500000000000), "Unexpected a0")
    require(diagonal == Q(40670769968163, 1000000000000), "Unexpected P(1)")

    radial = tuple(map(Q, data["radial_interval"]))
    angular = tuple(map(Q, data["target_angular_interval"]))
    radius = Q(data["kernel_support_radius"])
    cutoff = Q(data["spherical_cutoff"])
    require(radial == (Q(7, 5), Q(3, 2)), "Wrong radial interval")
    require(angular == (Q(-4, 5), Q(-11, 20)), "Wrong angular interval")
    require(radius == Q(7, 4) and cutoff == Q(7, 9), "Wrong support/cutoff")

    # d/dr [(r^2+s^2-1)/(2rs)] = (r^2-s^2+1)/(2r^2s).
    derivative_numerator_lower = radial[0]**2 - radial[1]**2 + 1
    require(derivative_numerator_lower > 0, "Radial monotonicity fails")
    require((2 * radial[1]**2 - 1) / (2 * radial[1]**2) == cutoff, "Cutoff identity fails")
    min_selected_distance_squared = 2 * radial[0]**2 * (1-angular[1])
    require(min_selected_distance_squared == Q(1519, 250) > radius**2, "Exterior gap fails")

    negative_p = [-c for c in power]
    negative_p_minus_one = list(negative_p)
    negative_p_minus_one[0] -= 1
    signs = [
        {"polynomial": "-P", **prove_strict_positive(negative_p, Q(-1), cutoff)},
        {"polynomial": "-P-1", **prove_strict_positive(negative_p_minus_one, *angular)},
    ]

    # For every real N, N*P(1)-a0*N^2 <= P(1)^2/(4*a0).
    real_maximum = diagonal**2 / (4 * harmonic[0])
    require(real_maximum < 1882, "Real-cardinality bound is insufficient")
    require(data["claimed_even_ordered_edge_bound"] == 1880, "Unexpected capacity")
    require(2 * (real_maximum.numerator // (2 * real_maximum.denominator)) == 1880,
            "Even ordered-edge rounding fails")
    center = diagonal / (2 * harmonic[0])
    integer_maximizers = (center.numerator // center.denominator,
                          center.numerator // center.denominator + 1)
    integer_maximum = max(n*diagonal-harmonic[0]*n*n for n in integer_maximizers)
    require(integer_maximum == Q(1881054173013741, 10**12), "Recorded integer maximum differs")

    return {
        "status": "INDEPENDENT_EXACT_MODERATE_GRAPH_CAPACITY_1880_PASS",
        "scope": "Local graph theorem only; does not by itself prove a global sphere-packing bound.",
        "input_sha256": hashlib.sha256(raw).hexdigest(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "arithmetic": "Python standard-library Fraction, exact rational polynomial division, primitive integer Sturm sequence",
        "archived_code_imported": False,
        "harmonic_degree": degree,
        "harmonic_coefficients": list(map(str, harmonic)),
        "independently_reconstructed_power_coefficients": list(map(str, power)),
        "closed_radial_interval": list(map(str, radial)),
        "closed_selected_angular_interval": list(map(str, angular)),
        "spherical_cutoff": str(cutoff),
        "radial_derivative_numerator_lower": str(derivative_numerator_lower),
        "whole_interval_sturm_proofs": signs,
        "a0": str(harmonic[0]),
        "P_at_1": str(diagonal),
        "real_cardinality_quadratic_upper_bound": str(real_maximum),
        "exact_gap_to_1882": str(1882-real_maximum),
        "integer_cardinality_maximum": str(integer_maximum),
        "even_ordered_edge_capacity": 1880,
        "minimum_selected_distance_squared": str(min_selected_distance_squared),
        "minimum_exterior_distance_gap_squared": str(min_selected_distance_squared-radius**2),
        "local_energy": "delta*(1880-ordered_selected_edge_count) >= 0 for every finite neighborhood and delta >= 0",
        "weight": "1/100000000",
        "density_charge": str(Q(1880, 8*10**8)),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    report = verify(args.certificate)
    text = json.dumps(report, indent=2) + "\n"
    if args.out:
        args.out.write_text(text)
    else:
        print(text, end="")
