"""Exact finite-code obstructions. No solver or floating arithmetic in proofs."""
from pathlib import Path
from fractions import Fraction
import argparse
import hashlib
import json
import math
import numpy as np

P = Path(__file__).resolve().parent
MUL = np.array([[0, 0, 0, 0], [0, 1, 2, 3],
                [0, 2, 3, 1], [0, 3, 1, 2]], dtype=np.uint8)


def shorten(G, column):
    G = G.copy()
    pivot = next(i for i in range(len(G)) if G[i, column])
    G[[0, pivot]] = G[[pivot, 0]]
    inverse = next(a for a in range(1, 4) if MUL[a, G[0, column]] == 1)
    G[0] = MUL[inverse, G[0]]
    for i in range(1, len(G)):
        G[i] ^= MUL[G[i, column], G[0]]
    assert not G[1:, column].any()
    return np.delete(G[1:], column, axis=1)


def verify(G, output):
    k, n = G.shape
    C = np.zeros((1, n), dtype=np.uint8)
    for row in G:
        C = np.concatenate([C ^ MUL[a, row] for a in range(4)])
    weights = np.count_nonzero(C, axis=1).astype(np.uint8)
    histogram = np.bincount(weights, minlength=n + 1)
    assert histogram[0] == 1  # Hence the linear encoding is injective.
    assert weights[1:].min() == 6
    forbidden = (weights > 0) & (weights < 8)
    spectrum = forbidden.astype(np.int64)
    degree = int(spectrum.sum())
    for bit in range(2 * k):
        blocks = spectrum.reshape(-1, 2, 1 << bit)
        a, b = blocks[:, 0, :].copy(), blocks[:, 1, :].copy()
        blocks[:, 0, :], blocks[:, 1, :] = a + b, a - b
    least = int(spectrum.min())
    assert int(spectrum[0]) == degree
    assert int(spectrum.sum()) == 0
    assert int(spectrum @ spectrum) == len(C) * degree
    assert (degree, least) == {17: (3330, -166), 18: (5271, -241)}[n]
    bound = Fraction(len(C) * -least, degree - least)
    integer_bound = bound.numerator // bound.denominator
    if n == 17:
        # Record delta_34 = 3^(33/2)/2^25; delta = M/3^(19/2).
        required = Fraction(3**26, 2**25)
        ratio = Fraction(integer_bound, 1) / required
    else:
        # New 36D delta = 2^24/(3^3 * 7^6); delta = M/3^10.
        required = Fraction(2**24 * 3**7, 7**6)
        ratio = Fraction(integer_bound, 1) / required
    assert integer_bound < required
    bad_words = C[forbidden]
    masks = np.unique(np.sum((bad_words != 0).astype(np.uint32)
                            * (1 << np.arange(n, dtype=np.uint32)), axis=1))
    assert len(masks) == {17: 1110, 18: 1757}[n]
    (output / f"supports{n}.txt").write_text(
        f"{n} {len(masks)}\n" + "\n".join(str(int(x)) for x in masks) + "\n")
    values, counts = np.unique(spectrum, return_counts=True)
    result = {
        "outer_length": n, "outer_dimension": k, "code_size": len(C),
        "packing_dimension": 2 * n,
        "weight_enumerator": histogram.tolist(), "forbidden_degree": degree,
        "least_character_eigenvalue": least,
        "independent_set_bound_exact": str(bound),
        "independent_set_bound_integer": integer_bound,
        "code_size_needed_to_beat_reference": required.numerator // required.denominator + 1,
        "density_ratio_upper_exact": str(ratio),
        "density_ratio_upper_decimal": float(ratio),
        "bad_support_count": len(masks),
        "character_eigenvalue_histogram": {str(int(v)): int(c) for v, c in zip(values, counts)},
        "weight_array_sha256": hashlib.sha256(weights.tobytes()).hexdigest(),
        "spectrum_sha256_le_i64": hashlib.sha256(spectrum.astype('<i8').tobytes()).hexdigest(),
        "scope": "Arbitrary distance-at-least-8 subsets of this specific outer code; full ternary parity fiber only.",
    }
    (output / f"spectral{n}.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in ["packing_dimension", "code_size",
          "forbidden_degree", "least_character_eigenvalue", "independent_set_bound_integer",
          "code_size_needed_to_beat_reference", "bad_support_count"]}), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=P / "work")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    G = np.loadtxt(P / "inputs/outer18_G11.txt", dtype=np.uint8)
    assert G.shape == (11, 18) and G.max() <= 3
    verify(shorten(G, 17), args.output_dir)
    verify(G, args.output_dir)


if __name__ == "__main__":
    main()
