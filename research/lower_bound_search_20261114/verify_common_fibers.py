"""Exact bounds on arbitrary common nonlinear fibers at squared minimum eight.

This covers the archived [20,13,6] outer code and all 1350 shortenings at
one, two or three coordinates. It does not constrain outer-word-dependent
fibers, different outer codes, a different minimum normalization, or metric changes.
"""
from pathlib import Path
from fractions import Fraction
import argparse
import hashlib
import itertools
import json
import math
import time
import numpy as np
from verify_codes import MUL

P = Path(__file__).resolve().parent
WEIGHTS = np.array([120, 5, 1], dtype=np.int64)
TARGETS = {
    17: Fraction(3**42, 2**45),
    18: Fraction(4 * 3**24, 7**6),
    19: Fraction(3**44, 5 * 2**43),
    20: Fraction(3**46, 2**46),
}


def outer_supports():
    G = np.loadtxt(P / 'inputs/outer20_G13.txt', dtype=np.uint8)
    assert G.shape == (13, 20)
    assert np.array_equal(G[:, :13], np.eye(13, dtype=np.uint8))
    C = np.zeros((1, 20), dtype=np.uint8)
    for row in G[:10]:
        C = np.concatenate([C ^ MUL[a, row] for a in range(4)])
    histogram = np.zeros(21, dtype=np.int64)
    supports = set()
    for block in range(64):
        offset = np.zeros(20, dtype=np.uint8)
        value = block
        for row in G[10:]:
            offset ^= MUL[value % 4, row]
            value //= 4
        words = C ^ offset
        weights = np.count_nonzero(words, axis=1)
        histogram += np.bincount(weights, minlength=21)
        short = words[(weights == 6) | (weights == 7)]
        masks = np.sum((short != 0).astype(np.uint32)
                       * (1 << np.arange(20, dtype=np.uint32)), axis=1)
        supports.update(map(int, masks))
    assert histogram.sum() == 4**13 and histogram[0] == 1
    assert not histogram[1:6].any()
    assert (int(histogram[6]), int(histogram[7])) == (3060, 9240)
    assert len(supports) == 4100
    return G, supports, histogram


def verified_group(supports):
    generators = json.loads((P / 'inputs/code20_support_automorphisms.json').read_text())['generators']
    for permutation in generators:
        assert sorted(permutation) == list(range(20))
        assert {sum(1 << permutation[i] for i in range(20) if support >> i & 1)
                for support in supports} == supports
    identity = tuple(range(20))
    group, pending = {identity}, [identity]
    for permutation in pending:
        for generator in generators:
            image = tuple(generator[permutation[i]] for i in range(20))
            if image not in group:
                group.add(image)
                pending.append(image)
    assert len(group) == 120
    return group


def column_rank(G, selected):
    A = G[:, selected].copy()
    rank = 0
    for column in range(len(selected)):
        choices = np.flatnonzero(A[rank:, column])
        if not len(choices):
            continue
        pivot = rank + int(choices[0])
        A[[rank, pivot]] = A[[pivot, rank]]
        inverse = next(a for a in (1, 2, 3) if MUL[a, A[rank, column]] == 1)
        A[rank] = MUL[inverse, A[rank]]
        for row in range(rank + 1, len(A)):
            A[row] ^= MUL[A[row, column], A[rank]]
        rank += 1
    return rank


def spectral_bound(n, supports):
    F = np.zeros((1 << n, 3), dtype=np.int64)
    for i in range(n):
        F[1 << i, 0] = 1
    for support in supports:
        weight = support.bit_count()
        assert weight in (6, 7)
        F[support, 1 if weight == 6 else 2] = 1
    for bit in range(n):
        blocks = F.reshape(-1, 2, 1 << bit, 3)
        a, b = blocks[:, 0].copy(), blocks[:, 1].copy()
        blocks[:, 0], blocks[:, 1] = a + 2*b, a - b
    degrees = F[0].copy()
    spectrum = F @ WEIGHTS
    degree, least = int(spectrum[0]), int(spectrum.min())
    # A character with nonzero-coordinate support J has multiplicity 2^|J|.
    multiplicity = np.array([1 << j.bit_count() for j in range(1 << n)], dtype=np.int64)
    assert int(multiplicity.sum()) == 3**n
    assert int(multiplicity @ spectrum) == 0
    assert int(multiplicity @ (spectrum*spectrum)) == 3**n * int(degrees @ (WEIGHTS*WEIGHTS))
    # Independent direct character sums at fixed, reproducible characters.
    rng = np.random.default_rng(20261114 + n)
    for character in [0, (1 << n)-1] + rng.integers(0, 1 << n, 30).tolist():
        direct = int(WEIGHTS[0]) * (2*n - 3*character.bit_count())
        for support in supports:
            k = (support & character).bit_count()
            w = support.bit_count()
            direct += int(WEIGHTS[1 if w == 6 else 2]) * (-1 if k % 2 else 1) * 2**(w-k)
        assert direct == int(spectrum[character])
    bound = Fraction(3**n * (-least), degree-least)
    ratio = bound / TARGETS[n]
    assert 0 < ratio < 1
    return {
        'packing_dimension': 2*n,
        'squared_minimum_target': 8,
        'short_support_count': len(supports),
        'edge_weights_for_ternary_weight1_outer_weight6_outer_weight7': WEIGHTS.tolist(),
        'weighted_degree': degree,
        'least_eigenvalue': least,
        'fiber_cardinality_upper_exact': str(bound),
        'fiber_cardinality_upper_floor': bound.numerator // bound.denominator,
        'fiber_cardinality_needed_to_beat_record': TARGETS[n].numerator // TARGETS[n].denominator + 1,
        'density_ratio_upper_exact': str(ratio),
        'density_ratio_upper_decimal': float(ratio),
        'spectrum_sha256_le_i64': hashlib.sha256(spectrum.astype('<i8').tobytes()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=P / 'work/common_fibers.json')
    args = parser.parse_args()
    start = time.time()
    G, supports, histogram = outer_supports()
    group = verified_group(supports)
    rows = []
    for count in (0, 1, 2, 3):
        unseen = set(itertools.combinations(range(20), count))
        covered = 0
        while unseen:
            selected = min(unseen)
            orbit = {tuple(sorted(permutation[i] for i in selected)) for permutation in group}
            assert orbit <= unseen
            unseen -= orbit
            covered += len(orbit)
            # The group need only preserve short supports. Verify the code size
            # separately for every shortening, even if the code groups differ.
            for element in orbit:
                assert column_rank(G, element) == count
            removed = sum(1 << i for i in selected)
            remaining = [i for i in range(20) if i not in selected]
            shortened = {sum(1 << j for j, i in enumerate(remaining) if support >> i & 1)
                         for support in supports if support & removed == 0}
            row = spectral_bound(20-count, shortened)
            row.update(removed_coordinates=list(selected), orbit_size=len(orbit), outer_dimension=13-count)
            rows.append(row)
            print(json.dumps(row), 'elapsed', round(time.time()-start, 2), flush=True)
        assert covered == math.comb(20, count)
    assert len(rows) == 18 and sum(row['orbit_size'] for row in rows) == 1351
    result = {'scope': __doc__, 'outer_weight_enumerator': histogram.tolist(),
              'verified_support_group_order': 120, 'outer_codes_covered': 1351,
              'shortenings_covered': 1350, 'orbit_results': rows}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print('PASS: arbitrary common nonlinear fibers fail the minimum-eight density targets for all 1351 outer codes.')


if __name__ == '__main__':
    main()
