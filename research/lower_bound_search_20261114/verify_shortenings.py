"""Exact rejection of all 1350 one-/two-/three-coordinate shortenings.

The supplied permutations are treated only as proposed symmetries. Their action
on every bad support is verified. No graph-canonicalization package is needed.
"""
from pathlib import Path
import argparse
import itertools
import json
import math
import subprocess
import tempfile
import numpy as np
from verify_codes import MUL

P = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=P / 'work/shortenings.json')
    args = parser.parse_args()
    G = np.loadtxt(P / 'inputs/outer20_G13.txt', dtype=np.uint8)
    assert G.shape == (13, 20) and np.array_equal(G[:, :13], np.eye(13, dtype=np.uint8))
    C = np.zeros((1, 20), dtype=np.uint8)
    for row in G[:10]:
        C = np.concatenate([C ^ MUL[a, row] for a in range(4)])
    histogram = np.zeros(21, dtype=np.int64)
    supports = []
    for index in range(64):
        offset = np.zeros(20, dtype=np.uint8)
        coefficient = index
        for row in G[10:]:
            offset ^= MUL[coefficient % 4, row]
            coefficient //= 4
        words = C ^ offset
        weights = np.count_nonzero(words, axis=1)
        histogram += np.bincount(weights, minlength=21)
        short = words[(weights == 6) | (weights == 7)]
        supports.extend(np.sum((short != 0).astype(np.uint32)
                               * (1 << np.arange(20, dtype=np.uint32)), axis=1).tolist())
    assert histogram.sum() == 4**13 and histogram[0] == 1 and not histogram[1:6].any()
    assert (int(histogram[6]), int(histogram[7])) == (3060, 9240)
    supports = set(map(int, supports))
    assert len(supports) == 4100
    generators = json.loads((P / 'inputs/code20_support_automorphisms.json').read_text())['generators']
    for permutation in generators:
        assert sorted(permutation) == list(range(20))
        images = {sum(1 << permutation[i] for i in range(20) if support >> i & 1)
                  for support in supports}
        assert images == supports
    # Independently compute the generated group's order.
    identity = tuple(range(20))
    group = {identity}
    pending = [identity]
    for permutation in pending:
        for generator in generators:
            image = tuple(generator[permutation[i]] for i in range(20))
            if image not in group:
                group.add(image)
                pending.append(image)
    assert len(group) == 120

    results = []
    with tempfile.TemporaryDirectory(prefix='shortening-replay-') as directory:
        W = Path(directory)
        binary = W / 'fiber_csp'
        subprocess.run(['g++', '-std=c++17', '-O3', str(P / 'fiber_csp.cpp'), '-o', str(binary)], check=True)
        for count in (1, 2, 3):
            unseen = set(itertools.combinations(range(20), count))
            covered = 0
            while unseen:
                selected = min(unseen)
                orbit = {tuple(sorted(permutation[i] for i in selected)) for permutation in group}
                assert orbit <= unseen
                unseen -= orbit
                covered += len(orbit)
                removed = sum(1 << i for i in selected)
                remaining = [i for i in range(20) if i not in selected]
                bad = sorted(sum(1 << j for j, i in enumerate(remaining) if support >> i & 1)
                             for support in supports if support & removed == 0)
                instance = W / 'instance.txt'
                instance.write_text(f'{len(remaining)} {len(bad)}\n' + '\n'.join(map(str, bad)) + '\n')
                result = subprocess.run([str(binary), str(instance), '1000000'],
                                        capture_output=True, text=True, check=True)
                assert result.stdout.startswith('UNSAT nodes '), result.stdout
                row = {'removed_coordinates': list(selected), 'orbit_size': len(orbit),
                       'packing_dimension': 2 * len(remaining), 'bad_supports': len(bad),
                       'result': result.stdout.strip()}
                results.append(row)
                print(json.dumps(row), flush=True)
            assert covered == math.comb(20, count)
    assert len(results) == 17 and sum(row['orbit_size'] for row in results) == 1350
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({'verified_support_group_order': len(group),
                                     'outer_weight_enumerator': histogram.tolist(),
                                     'shortenings_covered': 1350, 'orbit_results': results}, indent=2) + '\n')
    print('PASS: all 1350 shortenings rejected, using 17 verified permutation orbits.')


if __name__ == '__main__':
    main()
