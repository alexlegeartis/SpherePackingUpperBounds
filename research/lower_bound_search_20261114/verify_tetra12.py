"""Exact necessary-condition obstruction to products of ternary tetracodes.

This is a separate minimum-twelve target. It does not follow from the
minimum-eight common-fiber bounds. It does not exclude general self-orthogonal
ternary fibers or general higher-minimum constructions.
"""
from pathlib import Path
import argparse
import itertools
import json
import math
import subprocess
import tempfile
import numpy as np
from verify_common_fibers import outer_supports, verified_group

P = Path(__file__).resolve().parent


def write_instance(path, n, target, supports):
    six = sorted(s for s in supports if s.bit_count() == 6)
    seven = sorted(s for s in supports if s.bit_count() == 7)
    path.write_text(f'{n} {target} {len(six)} {len(seven)} 0\n'
                    + '\n'.join(map(str, six+seven)) + '\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=P / 'work/tetra12.json')
    args = parser.parse_args()
    _, supports, _ = outer_supports()
    group = verified_group(supports)
    rows = []
    with tempfile.TemporaryDirectory(prefix='tetra12-replay-') as directory:
        W = Path(directory)
        binary = W / 'search'
        subprocess.run(['g++', '-std=c++17', '-O3',
                        str(P / 'experiments/tetra12_search.cpp'), '-o', str(binary)], check=True)
        instance, output = W / 'instance.txt', W / 'result.txt'
        # Positive control, independently checked as two disjoint size-four
        # blocks. A negative control forbids every possible two-triple support.
        write_instance(instance, 8, 2, set())
        subprocess.run([str(binary), str(instance), str(output)], check=True, capture_output=True)
        lines = output.read_text().splitlines()
        assert lines[0].startswith('SAT ')
        blocks = list(map(int, lines[1:]))
        assert len(blocks) == 2 and all(b.bit_count() == 4 for b in blocks)
        assert not blocks[0] & blocks[1] and blocks[0] | blocks[1] == 255
        all_sixes = {sum(1 << i for i in selected) for selected in itertools.combinations(range(8), 6)}
        write_instance(instance, 8, 2, all_sixes)
        subprocess.run([str(binary), str(instance), str(output)], check=True, capture_output=True)
        assert output.read_text().startswith('UNSAT ')
        for count in (0, 1, 2, 3):
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
                shortened = {sum(1 << j for j, i in enumerate(remaining) if support >> i & 1)
                             for support in supports if support & removed == 0}
                n = len(remaining)
                target = 5 if n == 20 else 4
                write_instance(instance, n, target, shortened)
                subprocess.run([str(binary), str(instance), str(output)], check=True, capture_output=True)
                result = output.read_text().strip()
                assert result.startswith('UNSAT '), result
                row = {'packing_dimension': 2*n, 'removed_coordinates': list(selected),
                       'outer_codes_covered': len(orbit), 'required_tetracode_blocks': target,
                       'result': result}
                rows.append(row)
                print(json.dumps(row), flush=True)
            assert covered == math.comb(20, count)
    assert len(rows) == 18 and sum(row['outer_codes_covered'] for row in rows) == 1351
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({'scope': __doc__, 'outer_codes_covered': 1351,
                                     'orbit_results': rows}, indent=2) + '\n')
    print('PASS: all 1351 outer codes fail the necessary tetracode block-pair test.')


if __name__ == '__main__':
    main()
