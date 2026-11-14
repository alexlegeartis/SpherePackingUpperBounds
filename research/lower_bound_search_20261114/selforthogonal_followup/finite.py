"""Integer arithmetic shared by the follow-up searches and exact replays."""

from pathlib import Path
from fractions import Fraction
import itertools
import json
import sys
import numpy as np

P = Path(__file__).resolve().parent
ROOT = P.parent
sys.path.insert(0, str(ROOT))
from verify_codes import MUL, shorten

RECORD_SQUARED = {
    30: Fraction(1),
    32: Fraction(3**32, 2**48),
    34: Fraction(3**33, 2**50),
    36: Fraction(2**48, 3**6 * 7**12),
    38: Fraction(3**31, 25 * 2**38),
    40: Fraction(3**32, 2**40),
}


def ternary_code(G):
    C = np.zeros((1, G.shape[1]), dtype=np.int16)
    for row in G:
        C = np.concatenate([(C + a * row) % 3 for a in range(3)])
    assert np.count_nonzero(np.count_nonzero(C, axis=1) == 0) == 1
    return C


def support_masks(C):
    return np.sum(
        (C != 0).astype(np.uint32) * (1 << np.arange(C.shape[1], dtype=np.uint32)),
        axis=1,
    )


def outer_generator(n, removed=None):
    G = np.loadtxt(ROOT / "inputs/outer20_G13.txt", dtype=np.uint8)
    if removed is None:
        removed = list(range(n, 20))
    assert len(removed) == 20 - n
    for column in sorted(removed, reverse=True):
        G = shorten(G, column)
    assert G.shape == (n - 7, n)
    return G


def parent_supports():
    G = np.loadtxt(ROOT / "inputs/outer20_G13.txt", dtype=np.uint8)
    assert G.shape == (13, 20)
    assert np.array_equal(G[:, :13], np.eye(13, dtype=np.uint8))
    dual = np.c_[G[:, 13:].T, np.eye(7, dtype=np.uint8)]
    D = np.zeros((1, 20), dtype=np.uint8)
    for row in dual:
        D = np.concatenate([D ^ MUL[a, row] for a in range(4)])
    dual_minimum = int(np.count_nonzero(D[1:], axis=1).min())
    assert dual_minimum == 8  # All sets of at most seven columns are independent.
    C = np.zeros((1, 20), dtype=np.uint8)
    for row in G[:10]:
        C = np.concatenate([C ^ MUL[a, row] for a in range(4)])
    supports = {weight: set() for weight in (6, 7, 9)}
    histogram = np.zeros(21, dtype=np.int64)
    for block in range(64):
        offset = np.zeros(20, dtype=np.uint8)
        value = block
        for row in G[10:]:
            offset ^= MUL[value % 4, row]
            value //= 4
        words = C ^ offset
        weight = np.count_nonzero(words, axis=1)
        histogram += np.bincount(weight, minlength=21)
        for w in supports:
            supports[w].update(map(int, support_masks(words[weight == w])))
    assert histogram.sum() == 4**13 and histogram[0] == 1
    assert not histogram[1:6].any()
    assert [int(histogram[w]) for w in (6, 7, 9)] == [3060, 9240, 212400]
    assert [len(supports[w]) for w in (6, 7, 9)] == [1020, 3080, 65240]
    return supports, histogram, dual_minimum


def verified_group(supports):
    source = ROOT / "inputs/code20_support_automorphisms.json"
    generators = json.loads(source.read_text())["generators"]
    for permutation in generators:
        assert sorted(permutation) == list(range(20))
        for values in supports.values():
            A = np.array(sorted(values), dtype=np.uint32)
            B = np.zeros_like(A)
            for i in range(20):
                B |= ((A >> i) & 1) << permutation[i]
            assert np.array_equal(np.sort(B), A)
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


def shortening_orbits(group, count):
    unseen = set(itertools.combinations(range(20), count))
    while unseen:
        selected = min(unseen)
        orbit = {
            tuple(sorted(permutation[i] for i in selected)) for permutation in group
        }
        assert orbit <= unseen
        unseen -= orbit
        yield selected, len(orbit)


def shortened_supports(supports, selected):
    removed = sum(1 << i for i in selected)
    remaining = [i for i in range(20) if i not in selected]
    result = {}
    for weight, values in supports.items():
        A = np.array(sorted(values), dtype=np.uint32)
        A = A[A & removed == 0]
        B = np.zeros_like(A)
        for j, i in enumerate(remaining):
            B |= ((A >> i) & 1) << j
        result[weight] = set(map(int, B))
    return result


def write_instance(path, n, k, supports):
    path.write_text(
        f"{n} {k} "
        + " ".join(str(len(supports[w])) for w in (6, 7, 9))
        + "\n"
        + "\n".join(str(s) for w in (6, 7, 9) for s in sorted(supports[w]))
        + "\n"
    )


def audit_selforthogonal(G, outer_supports):
    assert np.all((0 <= G) & (G < 3))
    assert not (G @ G.T % 3).any()
    n, t = G.shape[1], len(G)
    C = ternary_code(G)
    weight = np.count_nonzero(C, axis=1)
    assert not (weight % 3).any()
    masks = support_masks(C)
    six, seven, nine = (outer_supports[w] for w in (6, 7, 9))
    contained = set(six)
    for support in seven:
        for i in range(n):
            if support >> i & 1:
                contained.add(support ^ (1 << i))
    bad6 = sum(int(s) in six for s in masks)
    bad9 = sum(int(s) in contained - six or int(s) in nine for s in masks)
    minimum = 6 if bad6 else 9 if bad9 else 12
    delta12_squared = Fraction(2) ** (2 * n - 28) * Fraction(3) ** (2 * t - n)
    actual_squared = delta12_squared * Fraction(minimum, 12) ** (2 * n)
    return {
        "dimension": 2 * n,
        "outer_dimension": n - 7,
        "ternary_dimension": t,
        "attained_minimum_squared": minimum,
        "norm6_bad_trit_words": bad6,
        "norm9_bad_trit_words": bad9,
        "search_energy": (10 * bad6 + bad9) // 2,
        "conditional_minimum12_density_squared": str(delta12_squared),
        "actual_density_squared": str(actual_squared),
        "actual_density_ratio_to_pinned_record_squared": str(
            actual_squared / RECORD_SQUARED[2 * n]
        ),
        "beats_pinned_record": actual_squared > RECORD_SQUARED[2 * n],
    }


def geometry_costs(n):
    """Exact minimum cost for each trit support; used for floating-free search."""
    assert n in (16, 17, 18)
    G = outer_generator(n)
    C = np.zeros((1, n), dtype=np.uint8)
    for row in G:
        C = np.concatenate([C ^ MUL[a, row] for a in range(4)])
    masks = support_masks(C)
    weight = np.array([s.bit_count() for s in range(1 << n)], dtype=np.int16)
    distance = np.full(1 << n, 12, dtype=np.int16)
    distance[masks] = np.minimum(weight[masks], 12)
    for i in range(n):
        block = distance.reshape(-1, 2, 1 << i)
        a, b = block[:, 0].copy(), block[:, 1].copy()
        block[:, 0], block[:, 1] = np.minimum(a, b + 2), np.minimum(b, a + 4)
    return distance
