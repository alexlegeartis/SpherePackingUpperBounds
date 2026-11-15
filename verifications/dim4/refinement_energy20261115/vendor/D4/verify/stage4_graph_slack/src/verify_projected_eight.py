#!/usr/bin/env python3
"""Replay the exact projected S^2 code certificate.

The certificate is a nonnegative combination of normalized Legendre kernels.
This checker reconstructs the polynomial with rational arithmetic, proves its
sign on the complete closed interval by Bernstein subdivision, and derives
the integer code bound from the exact objective.  No samples or floats enter
any decision.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import json
import math
import sys

sys.set_int_max_str_digits(0)

EXPECTED_CUTOFF = Q(223, 735)
EXPECTED_CONSTANT = Q(99999, 100000)
EXPECTED_DENOMINATOR = 10_000_000_000
MAX_DEPTH = 20
MAX_NODES = 10_000


def parse_exact(value):
    if type(value) not in (str, int):
        raise ValueError(f"expected an exact rational string or integer, got {value!r}")
    return Q(value)


def add(a, b):
    result = [Q(0)] * max(len(a), len(b))
    for i, value in enumerate(a):
        result[i] += value
    for i, value in enumerate(b):
        result[i] += value
    return result


def scale(a, factor):
    return [factor * value for value in a]


def legendre_power(degree):
    """Return normalized P_degree in the monomial basis."""
    polynomials = [[Q(1)], [Q(0), Q(1)]]
    for ell in range(1, degree):
        # P_(ell+1) = ((2 ell + 1) t P_ell - ell P_(ell-1))/(ell+1).
        next_poly = [Q(0)] * (ell + 2)
        for i, value in enumerate(polynomials[-1]):
            next_poly[i + 1] += Q(2 * ell + 1, ell + 1) * value
        for i, value in enumerate(polynomials[-2]):
            next_poly[i] -= Q(ell, ell + 1) * value
        polynomials.append(next_poly)
    return polynomials[degree]


def certificate_polynomial(coefficients):
    result = [Q(0)]
    for degree, coefficient in enumerate(coefficients):
        result = add(result, scale(legendre_power(degree), coefficient))
    return result


def compose_affine(polynomial, slope, intercept):
    """Return polynomial(u) = polynomial(slope*u + intercept)."""
    result = [Q(0)] * len(polynomial)
    for degree, coefficient in enumerate(polynomial):
        for power in range(degree + 1):
            result[power] += (coefficient * math.comb(degree, power)
                              * slope ** power * intercept ** (degree - power))
    return result


def bernstein_coefficients(polynomial, left, width, degree):
    """Bernstein coefficients on [left, left + width]."""
    # First compose with x = left + width*z, then convert powers of z to
    # degree-degree Bernstein coefficients exactly.
    powers = [sum(polynomial[j] * math.comb(j, power)
                  * left ** (j - power) * width ** power
                  for j in range(power, degree + 1))
              for power in range(degree + 1)]
    return [sum(powers[power] * Q(math.comb(index, power), math.comb(degree, power))
                for power in range(index + 1))
            for index in range(degree + 1)]


def prove_interval(polynomial):
    """Prove polynomial >= 0 on [0,1] by exact dyadic Bernstein cells."""
    degree = len(polynomial) - 1
    stack = [(0, 0)]
    leaves = []
    nodes = 0
    while stack:
        numerator, depth = stack.pop()
        nodes += 1
        if nodes > MAX_NODES:
            raise ValueError("Bernstein subdivision exceeded node budget")
        denominator = 1 << depth
        coefficients = bernstein_coefficients(
            polynomial, Q(numerator, denominator), Q(1, denominator), degree)
        if all(value >= 0 for value in coefficients):
            leaves.append((numerator, depth))
            continue
        if depth >= MAX_DEPTH:
            raise ValueError("Bernstein subdivision did not certify a leaf")
        stack.append((2 * numerator + 1, depth + 1))
        stack.append((2 * numerator, depth + 1))
    return nodes, leaves


def verify(path):
    data = json.loads(Path(path).read_text())
    if data.get("cutoff") != str(EXPECTED_CUTOFF):
        raise ValueError("unexpected projected-code cutoff")
    if data.get("coefficient_denominator") != EXPECTED_DENOMINATOR:
        raise ValueError("unexpected coefficient denominator")
    coefficients = [parse_exact(value) for value in data["coefficients"]]
    if len(coefficients) != 28:
        raise ValueError("projected certificate must have degree 27")
    if coefficients[0] != EXPECTED_CONSTANT or coefficients[0] <= 0:
        raise ValueError("invalid constant Legendre coefficient")
    if any(value < 0 for value in coefficients[1:]):
        raise ValueError("Legendre coefficients must be nonnegative")

    polynomial = certificate_polynomial(coefficients)
    transformed = scale(compose_affine(polynomial, Q(958, 735), Q(-1)), Q(-1))
    nodes, leaves = prove_interval(transformed)
    endpoint = sum(transformed)
    if endpoint < 0:
        raise ValueError("projected polynomial sign failed at the interval endpoint")

    objective = sum(coefficients) / coefficients[0]
    if not objective < 9:
        raise ValueError("projected certificate objective is not strictly below nine")
    count_bound = objective.numerator // objective.denominator
    if count_bound != 8:
        raise ValueError(f"projected certificate gives unexpected count bound {count_bound}")

    return {
        "status": "EXACT_PROJECTED_EIGHT_PASS",
        "cutoff": str(EXPECTED_CUTOFF),
        "degree": len(coefficients) - 1,
        "constant_coefficient": str(coefficients[0]),
        "objective": str(objective),
        "objective_strictly_below_nine": True,
        "count_bound": count_bound,
        "bernstein_nodes": nodes,
        "bernstein_leaves": len(leaves),
        "bernstein_depth": max(depth for _, depth in leaves),
        "endpoint_check": str(endpoint),
        "sign_interval": "[-1,223/735]",
        "sign_method": "exact rational Bernstein subdivision",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path)
    args = parser.parse_args()
    try:
        result = verify(args.certificate)
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
        print(json.dumps({"status": "FAIL", "error": str(exc)}))
        return 1
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
