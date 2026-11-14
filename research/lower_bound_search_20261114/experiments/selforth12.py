"""Search or audit ternary self-orthogonal fibers for a minimum-twelve target.

Search is heuristic. Auditing a returned code and its attained minimum uses
only exact finite enumeration. A nonzero violation count is not a packing
with the conditional minimum-twelve density.
"""
from pathlib import Path
from fractions import Fraction
import argparse
import hashlib
import json
import subprocess
import sys
import numpy as np

P = Path(__file__).resolve().parent
sys.path.insert(0, str(P.parent))
from verify_codes import MUL, shorten


def support_masks(words):
    n = words.shape[1]
    return np.sum((words != 0).astype(np.uint32)
                  * (1 << np.arange(n, dtype=np.uint32)), axis=1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dimension', type=int, choices=(34, 36), required=True)
    parser.add_argument('--seconds', type=float, default=120)
    parser.add_argument('--steps', type=int, default=2000000)
    parser.add_argument('--output-prefix', type=Path)
    parser.add_argument('--audit-best', type=Path)
    args = parser.parse_args()
    n = args.dimension // 2
    G = np.loadtxt(P.parent / 'inputs/outer18_G11.txt', dtype=np.uint8)
    if n == 17:
        G = shorten(G, 17)
    C = np.zeros((1, n), dtype=np.uint8)
    for row in G:
        C = np.concatenate([C ^ MUL[a, row] for a in range(4)])
    weights = np.count_nonzero(C, axis=1)
    assert len(C) == 4**len(G) and weights[1:].min() == 6
    supports = [set(map(int, support_masks(C[weights == weight]))) for weight in (6, 7, 9)]
    six, seven, nine = supports
    bad_six = set(six)
    for support in seven:
        for i in range(n):
            if support >> i & 1:
                bad_six.add(support ^ (1 << i))
    if args.output_prefix is None:
        args.output_prefix = P.parent / f'work/selforth12_{n}'
    prefix = args.output_prefix.resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    if args.audit_best is None:
        initial = np.zeros((8, n), dtype=np.int16)
        for j in range(4):
            initial[2*j, 4*j:4*j+4] = [1, 1, 1, 0]
            initial[2*j+1, 4*j:4*j+4] = [0, 1, 2, 1]
        instance = Path(str(prefix) + '_input.txt')
        instance.write_text(f'{n} 8 ' + ' '.join(str(len(s)) for s in supports) + '\n'
                            + '\n'.join(str(s) for group in supports for s in sorted(group)) + '\n'
                            + '\n'.join(' '.join(map(str, row)) for row in initial) + '\n')
        binary = Path(str(prefix) + '_search')
        subprocess.run(['g++', '-std=c++17', '-O3', str(P / 'selforth12_search.cpp'),
                        '-o', str(binary)], check=True)
        best_file = Path(str(prefix) + '_best.txt')
        with Path(str(prefix) + '_search.log').open('w') as log:
            subprocess.run([str(binary), str(instance), str(best_file),
                            str(args.seconds), str(args.steps)], stdout=log, stderr=subprocess.STDOUT, check=True)
    else:
        best_file = args.audit_best
    lines = best_file.read_text().splitlines()
    header = list(map(int, lines[0].split()))
    assert header[:2] == [n, 8]
    H = np.array([list(map(int, line.split())) for line in lines[1:]], dtype=np.int16)
    assert H.shape == (8, n) and np.all((0 <= H) & (H < 3))
    assert not (H @ H.T % 3).any()
    T = np.zeros((1, n), dtype=np.int16)
    for row in H:
        T = np.concatenate([(T + a*row) % 3 for a in range(3)])
    tw = np.count_nonzero(T, axis=1)
    assert len(T) == 3**8 and np.count_nonzero(tw == 0) == 1
    assert not (tw % 3).any()
    tm = support_masks(T)
    norm6_count = sum(int(s) in six for s in tm)
    norm9_count = sum(int(s) in bad_six-six or int(s) in nine for s in tm)
    assert 10*norm6_count + norm9_count == 2*header[2]
    attained_minimum = 6 if norm6_count else 9 if norm9_count else 12
    # Self-orthogonality forces all squared norms to be multiples of three.
    # Weight-three trit differences are safe at twelve, and weights >= twelve
    # are automatically safe. Weight six can only give norm six or nine via
    # a weight-six match or containment in a weight-seven support; weight nine
    # can only give norm nine via a matching weight-nine outer support.
    # The period itself supplies norm twelve, so this minimum is attained.
    delta12_squared = Fraction(2)**(4*len(G)-2*n) * Fraction(3)**(16-n)
    actual_squared = delta12_squared * Fraction(attained_minimum, 12)**(2*n)
    record_squared = Fraction(3**33, 2**50) if n == 17 else Fraction(2**48, 3**6*7**12)
    result = {
        'packing_dimension': 2*n, 'outer_dimension': len(G), 'ternary_dimension': 8,
        'target_minimum_squared': 12, 'attained_minimum_squared': attained_minimum,
        'norm6_bad_trit_words': norm6_count, 'norm9_bad_trit_words': norm9_count,
        'conditional_center_density_at_minimum12_squared': str(delta12_squared),
        'actual_center_density_squared': str(actual_squared),
        'actual_density_ratio_to_pinned_record_squared': str(actual_squared / record_squared),
        'is_new_record_against_pinned_value': actual_squared > record_squared,
        'ternary_generator_sha256_u8': hashlib.sha256(H.astype(np.uint8).tobytes()).hexdigest(),
    }
    Path(str(prefix) + '_audit.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
