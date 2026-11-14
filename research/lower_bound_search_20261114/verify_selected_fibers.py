"""Exact bounds allowing outer-word-dependent selected-coordinate fibers.

Fix the archived [17,10,6] outer code. Fix at most four trit coordinates
and a parity check on all remaining trits. Allow arbitrary unions of those
ternary-fiber cosets, with arbitrary dependence on every outer word.
The exact conflict graph has group F2^20 x F3^k. Its eigenvalues are obtained
without materializing its millions of vertices or its adjacency matrix.
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
from verify_codes import shorten, MUL

P = Path(__file__).resolve().parent


def walsh(values):
    result = values.copy()
    for bit in range(result.size.bit_length() - 1):
        z = result.reshape(-1, 2, 1 << bit)
        a, b = z[:, 0, :].copy(), z[:, 1, :].copy()
        z[:, 0, :], z[:, 1, :] = a + b, a - b
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=P / 'work/selected_fibers.json')
    parser.add_argument('--max-selected', type=int, choices=(3, 4), default=3)
    args = parser.parse_args()
    G = shorten(np.loadtxt(P / 'inputs/outer18_G11.txt', dtype=np.uint8), 17)
    C = np.zeros((1, 17), dtype=np.uint8)
    for row in G:
        C = np.concatenate([C ^ MUL[a, row] for a in range(4)])
    weight = np.count_nonzero(C, axis=1)
    forbidden = (weight > 0) & (weight < 8)
    assert forbidden.sum() == 3330 and weight[1:].min() == 6
    base = walsh(forbidden.astype(np.int32))
    assert base.min() == -166
    N = len(C)
    indices = np.arange(N, dtype=np.int32)
    shifts = []
    for coordinate in range(17):
        low = high = 0
        for row in range(10):
            for bit, coefficient in enumerate((1, 2)):
                value = int(MUL[coefficient, G[row, coordinate]])
                low |= (value & 1) << (2 * row + bit)
                high |= ((value >> 1) & 1) << (2 * row + bit)
        shifts.append((low, high))

    def extend(F, coordinate, nonzero_character):
        low, high = shifts[coordinate]
        S = F + F[indices ^ low] + F[indices ^ high] + F[indices ^ low ^ high]
        if nonzero_character:
            assert not np.any(S % 2)
            return -F + S // 2
        assert not np.any(S % 4)
        return 2 * F - S // 4

    # Independently compare the shift formula with direct Walsh transforms.
    # Checking the encoding basis verifies every coordinate's character labels.
    for coordinate, (low, high) in enumerate(shifts):
        for bit in range(20):
            value = int(C[1 << bit, coordinate])
            assert (value & 1) == ((low >> bit) & 1)
            assert ((value >> 1) & 1) == ((high >> bit) & 1)
    for nonzero in (False, True):
        direct = forbidden.astype(np.int32)
        derived = base
        for coordinate in (0, 1, 2):
            direct *= np.where(C[:, coordinate] != 0, -1 if nonzero else 2, 1)
            derived = extend(derived, coordinate, nonzero)
        assert np.array_equal(walsh(direct), derived)
    del C

    rows = []
    start = time.time()

    def record(selected, spectra):
        degree = int(spectra[0][0])
        least = min(int(F.min()) for F in spectra)
        k = len(selected)
        bound = Fraction(N * 3**k * -least, degree - least)
        target = Fraction(3**(26 + k), 2**25)
        ratio = bound / target
        rows.append({'selected_coordinates': list(selected), 'degree': degree,
                     'least_eigenvalue': least, 'states_upper_bound': str(bound),
                     'states_needed_for_record': target.numerator // target.denominator + 1,
                     'density_ratio_upper_exact': str(ratio),
                     'density_ratio_upper_decimal': float(ratio)})
        assert ratio < 1, f'Unexcluded family at {selected}: {ratio}'

    for a in range(17):
        first = [extend(base, a, nonzero) for nonzero in (False, True)]
        record((a,), first)
        for b in range(a + 1, 17):
            second = [extend(F, b, nonzero) for F in first for nonzero in (False, True)]
            record((a, b), second)
            for c in range(b + 1, 17):
                third = [extend(F, c, nonzero) for F in second for nonzero in (False, True)]
                record((a, b, c), third)
                if args.max_selected == 4:
                    for d in range(c + 1, 17):
                        fourth = [extend(F, d, nonzero) for F in third for nonzero in (False, True)]
                        record((a, b, c, d), fourth)
        print(f'Completed first coordinate {a}; checked {len(rows)} selections; '
              f'elapsed {time.time() - start:.1f}s', flush=True)

    expected = sum(math.comb(17, k) for k in range(1, args.max_selected+1))
    assert len(rows) == expected
    summary = {}
    for k in range(1, args.max_selected+1):
        subset = [r for r in rows if len(r['selected_coordinates']) == k]
        summary[str(k)] = max(subset, key=lambda r: Fraction(r['density_ratio_upper_exact']))
    result = {'scope': __doc__, 'selections_checked': len(rows),
              'base_spectrum_sha256_le_i32': hashlib.sha256(base.astype('<i4').tobytes()).hexdigest(),
              'worst_density_bound_by_selected_count': summary, 'results': rows}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(summary, indent=2), flush=True)
    print(f'PASS: all {expected} selected-coordinate fiber families are below the pinned 34D record.')


if __name__ == '__main__':
    main()
