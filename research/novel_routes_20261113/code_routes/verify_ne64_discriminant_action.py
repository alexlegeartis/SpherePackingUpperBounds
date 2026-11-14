#!/usr/bin/env python3
"""Exact finite-group audit for a proposed Ne64 Eisenstein-unit argument.

Uses only NumPy integer arithmetic; no tolerance-based equality is used.
Run from any working directory.  Cached full groups are optional and are not
trusted by this checker: both factors are generated again from the source.
"""

from pathlib import Path
import argparse
import hashlib
import json
import re

import numpy as np


HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "dual64"


def parse_generators(path):
    text = path.read_text()
    section = text.split('<a NAME="SUBGROUP_GENERATORS">')[1]
    section = section.split('<a NAME="PROPERTIES">')[0]
    numbers = []
    for line in section.splitlines():
        line = re.sub('<[^>]*>', '', line).strip()
        if re.fullmatch(r'[0-9+\-\s]+', line or 'x'):
            numbers.extend(map(int, line.split()))
    assert numbers[0] == 6
    p = 1
    answer = []
    for _ in range(6):
        assert numbers[p:p+2] == [64, 64]
        p += 2
        answer.append(np.array(numbers[p:p+4096], dtype=np.int64).reshape(64, 64))
        p += 4096
    assert p == len(numbers)
    return answer


def product(a, b):
    # Explicitly exclude overflow in NumPy's exact integer multiplication.
    bound = int(a.shape[1]) * int(np.max(np.abs(a))) * int(np.max(np.abs(b)))
    assert bound < 2**63
    return a @ b


def generated_group(generators):
    identity = np.eye(64, dtype=np.int64)
    group = [identity]
    seen = {identity.tobytes()}
    cursor = 0
    while cursor < len(group):
        a = group[cursor]
        cursor += 1
        for g in generators:
            b = product(a, g)
            key = b.tobytes()
            if key not in seen:
                seen.add(key)
                group.append(b)
        assert len(group) <= 5000
    return group


def rank_mod3(matrix):
    a = matrix.copy() % 3
    r = 0
    for c in range(a.shape[1]):
        pivot = next((i for i in range(r, a.shape[0]) if a[i, c]), None)
        if pivot is None:
            continue
        a[[r, pivot]] = a[[pivot, r]]
        a[r] = a[r] * int(a[r, c]) % 3
        for i in range(a.shape[0]):
            if i != r and a[i, c]:
                a[i] = (a[i] - a[i, c] * a[r]) % 3
        r += 1
    return r


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--save-groups", action="store_true",
                        help="also regenerate the ignored full-group NumPy caches")
    args = parser.parse_args()
    source = np.load(DATA / "source.npz")
    G, H = source["G"], source["H"]
    identity = np.eye(64, dtype=np.int64)
    assert np.array_equal(product(G, H), 3 * identity)
    assert rank_mod3(H) == 32
    generators = parse_generators(HERE / "CQ64a.html")
    assert all(np.array_equal(product(product(g, G), g.T), G)
               for g in generators[:5])
    assert not np.array_equal(product(product(generators[5], G), generators[5].T), G)
    assert all(np.array_equal(product(a, b), product(b, a))
               for a in generators[:2] for b in generators[2:5])

    # Coordinates are rows.  A dual vector has coordinates z G^{-1}, and
    # g fixes every discriminant class iff G^{-1}(g-I) is integral.
    def fixes_discriminant(g):
        return bool(np.all(product(H, g - identity) % 3 == 0))

    large = generated_group(generators[:2])
    small = generated_group(generators[2:5])
    assert len(large) == 4896
    assert len(small) == 120
    if args.save_groups:
        np.savez_compressed(HERE / "sl2_17_group.npz", group=np.array(large, dtype=np.int64))
        np.savez_compressed(HERE / "sl2_5_group.npz", group=np.array(small, dtype=np.int64))
    kernels = [sum(fixes_discriminant(g) for g in group)
               for group in (large, small)]
    assert kernels == [1, 1]
    assert not fixes_discriminant(-identity)
    intersection = {g.tobytes() for g in large}.intersection(
        g.tobytes() for g in small)
    assert intersection == {identity.tobytes(), (-identity).tobytes()}

    tau = [g for g in small
           if np.array_equal(product(g, g) + g + identity, 0 * identity)]
    assert len(tau) == 20
    fixed_dimensions = sorted({32 - rank_mod3(product(H, g - identity))
                               for g in tau})
    assert fixed_dimensions == [16]
    assert all(not fixes_discriminant(g) for g in tau)

    # Only the identity and the order-three elements of this factor fix
    # any nonzero discriminant classes.  This is also checked directly.
    other_fixed_dimensions = sorted({32 - rank_mod3(product(H, g - identity))
                                     for g in small
                                     if not np.array_equal(g, identity)
                                     and not np.array_equal(
                                         product(g, g) + g + identity,
                                         0 * identity)})
    assert other_fixed_dimensions == [0]

    result = {
        "source_url": "https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/CQ64a.html",
        "source_sha256": hashlib.sha256((HERE / "CQ64a.html").read_bytes()).hexdigest(),
        "ne64_gram_source_sha256": hashlib.sha256((DATA / "CQ64b.html").read_bytes()).hexdigest(),
        "first_five_generators_preserve_G": True,
        "sixth_generator_preserves_G": False,
        "factor_orders": [4896, 120],
        "factor_intersection_order": 2,
        "central_product_order": 293760,
        "factor_discriminant_kernel_orders": kernels,
        "full_known_subgroup_discriminant_kernel_order": 1,
        "SL2_5_fixed_point_free_order3_elements": len(tau),
        "their_discriminant_fixed_dimensions": fixed_dimensions,
        "other_nonidentity_SL2_5_discriminant_fixed_dimensions": other_fixed_dimensions,
        "scope": "The published subgroup only; no assertion that it is the full automorphism group.",
    }
    (HERE / "ne64_discriminant_action_verified.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
