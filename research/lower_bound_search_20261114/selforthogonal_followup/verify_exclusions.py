"""Replay new finite exclusions. No SAT solver or floating arithmetic required."""

from pathlib import Path
import argparse
import itertools
import json
import math
import re
import subprocess
import sys
import tempfile
import numpy as np
from finite import (
    parent_supports,
    verified_group,
    shortening_orbits,
    shortened_supports,
    write_instance,
    ternary_code,
    support_masks,
)

P = Path(__file__).resolve().parent


def saved_factorizations():
    result = []
    cases = [
        (
            17,
            P.parent / "results/selforth12_17_best.txt",
            "factorization17.json",
            "original34",
        ),
        (
            18,
            P.parent / "results/selforth12_18_best.txt",
            "factorization18.json",
            "original36",
        ),
        (18, P / "results/neighbor18_best.txt", "factorization18_new.json", "new36"),
    ]
    for n, candidate, witness_file, label in cases:
        G = np.loadtxt(candidate, skiprows=1, dtype=np.int16)
        assert G.shape == (8, n) and not (G @ G.T % 3).any()
        T = ternary_code(G)
        witness = json.loads((P / "results" / witness_file).read_text())
        groups = witness["groups"]
        assert sorted(i for group in groups for i in group) == list(range(n))
        sizes = []
        for group in groups:
            C = np.unique(T[:, group], axis=0)
            sizes.append(len(C))
            weights = np.count_nonzero(C, axis=1)
            masks = support_masks(C)
            three = set(map(int, masks[weights == 3]))
            if len(group) == 4:
                assert len(C) == 9 and three == {15 ^ (1 << i) for i in range(4)}
            else:
                assert len(group) == 9 and len(C) == 81 and len(three) == 3
                assert all(not a & b for a, b in itertools.combinations(three, 2))
                assert sum(three) == 511
                expected = {a | b for a, b in itertools.combinations(three, 2)}
                a, b, c = sorted(three)
                for i in range(9):
                    if not a >> i & 1:
                        continue
                    for j in range(9):
                        if not b >> j & 1:
                            continue
                        for k in range(9):
                            if c >> k & 1:
                                expected.add(511 ^ (1 << i) ^ (1 << j) ^ (1 << k))
                assert len(expected) == 30
                assert set(map(int, masks[weights == 6])) == expected
                assert set(map(int, masks[weights == 9])) == {511}
        assert math.prod(sizes) == len(T) == 3**8
        result.append(
            {
                "dimension": 2 * n,
                "candidate": label,
                "factor_lengths": list(map(len, groups)),
                "factor_code_sizes": sizes,
                "direct_product_verified": True,
            }
        )
    return result


def golay_obstruction(supports):
    G = np.loadtxt(P / "inputs/golay12_G6.txt", dtype=np.int16)
    assert G.shape == (6, 12) and not (G @ G.T % 3).any()
    C = ternary_code(G)
    weights = np.count_nonzero(C, axis=1)
    assert weights[1:].min() == 6
    blocks = set(map(int, support_masks(C[weights == 6])))
    assert len(blocks) == 132
    for selected in itertools.combinations(range(12), 7):
        mask = sum(1 << i for i in selected)
        assert sum(s & mask == s for s in blocks) == 1
    contains = np.zeros(1 << 20, dtype=np.uint8)
    contains[list(supports[7])] = 1
    for i in range(20):
        block = contains.reshape(-1, 2, 1 << i)
        block[:, 1] |= block[:, 0]
    twelve = np.array(
        [s for s in range(1 << 20) if s.bit_count() == 12], dtype=np.uint32
    )
    assert len(twelve) == 125970 and np.all(contains[twelve])
    return {
        "golay_weight6_supports": 132,
        "golay_7_subsets_checked": 792,
        "parent_12_subsets_checked": len(twelve),
        "conclusion": "Every placement of a ternary Golay12 direct factor has a norm-nine violation.",
    }


def controls(directory, components, products):
    empty = {6: set(), 7: set(), 9: set()}
    instance, output = directory / "control.txt", directory / "control.out"
    write_instance(instance, 9, 4, empty)
    subprocess.run(
        [str(components), str(instance), str(output), "component"],
        check=True,
        capture_output=True,
    )
    lines = output.read_text().splitlines()
    assert lines[0].startswith("SAT ")
    G = np.array([list(map(int, line.split())) for line in lines[1:]], dtype=np.int16)
    assert G.shape == (4, 9) and not (G @ G.T % 3).any()
    T = ternary_code(G)
    assert np.bincount(np.count_nonzero(T, axis=1), minlength=10).tolist() == [
        1,
        0,
        0,
        6,
        0,
        0,
        66,
        0,
        0,
        8,
    ]
    write_instance(instance, 9, 4, {6: set(), 7: set(), 9: {511}})
    subprocess.run(
        [str(components), str(instance), str(output), "component"],
        check=True,
        capture_output=True,
    )
    assert output.read_text().startswith("UNSAT ")
    for n in (15, 16):
        write_instance(instance, n, 0, empty)
        subprocess.run(
            [str(products), str(instance), str(output)], check=True, capture_output=True
        )
        lines = output.read_text().splitlines()
        assert lines[0].startswith("SAT ")
        blocks = [int(line) for line in lines[1:] if not line.startswith("REPEAT")]
        assert len(blocks) == n // 4 and all(s.bit_count() == 4 for s in blocks)
        assert all(not a & b for a, b in itertools.combinations(blocks, 2))
        if n == 15:
            repeat = next(
                int(line.split()[1]) for line in lines[1:] if line.startswith("REPEAT")
            )
            assert repeat.bit_count() == 3 and repeat | sum(blocks) == (1 << n) - 1
        all_six = {
            sum(1 << i for i in selected)
            for selected in itertools.combinations(range(n), 6)
        }
        write_instance(instance, n, 0, {6: all_six, 7: set(), 9: set()})
        subprocess.run(
            [str(products), str(instance), str(output)], check=True, capture_output=True
        )
        assert output.read_text().startswith("UNSAT ")


def main():
    assert sys.flags.optimize == 0, "Do not disable assertions"
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=P / "work/exclusions.json")
    args = parser.parse_args()
    supports, histogram, dual_minimum = parent_supports()
    group = verified_group(supports)
    golay = golay_obstruction(supports)
    factorizations = saved_factorizations()
    rows = []
    with tempfile.TemporaryDirectory(prefix="selforth-followup-") as directory:
        W = Path(directory)
        components, products = W / "components", W / "products"
        for name, binary in [
            ("component_search.cpp", components),
            ("product_search.cpp", products),
        ]:
            subprocess.run(
                ["g++", "-std=c++17", "-O3", str(P / "src" / name), "-o", str(binary)],
                check=True,
            )
        controls(W, components, products)
        for count in (2, 3, 4, 5):
            covered = 0
            cases = 0
            for selected, orbit_size in shortening_orbits(group, count):
                current = shortened_supports(supports, selected)
                n = 20 - count
                instance, output = W / "instance.txt", W / "result.txt"
                write_instance(instance, n, 0, current)
                binary = components if count in (2, 3) else products
                subprocess.run(
                    [str(binary), str(instance), str(output)],
                    check=True,
                    capture_output=True,
                )
                result = output.read_text().strip()
                assert result.startswith("UNSAT "), result
                if count in (2, 3):
                    tested = int(re.search(r"tested_components (\d+)", result)[1])
                    assert tested == math.comb(n, 9) * 280
                    if n == 18:
                        assert "valid_components 0 " in result
                row = {
                    "dimension": 2 * n,
                    "removed_coordinates": list(selected),
                    "outer_codes_covered": orbit_size,
                    "fiber_family": (
                        "T9 direct factor"
                        if n == 18
                        else (
                            "T9 + two tetracodes"
                            if n == 17
                            else (
                                "four tetracodes"
                                if n == 16
                                else "three tetracodes + repetition3"
                            )
                        )
                    ),
                    "result": result,
                }
                rows.append(row)
                covered += orbit_size
                cases += 1
            assert covered == math.comb(20, count)
            assert cases == {2: 4, 3: 12, 4: 53, 5: 137}[count]
            print(
                f"PASS shortening {count}: {cases} orbits, {covered} outer codes.",
                flush=True,
            )
    assert len(rows) == 206 and sum(row["outer_codes_covered"] for row in rows) == 21679
    result = {
        "scope": "Only the explicitly listed minimum-twelve fiber families.",
        "outer_weight_enumerator": histogram.tolist(),
        "outer_dual_minimum": dual_minimum,
        "verified_support_group_order": len(group),
        "golay_obstruction": golay,
        "saved_seed_factorizations": factorizations,
        "outer_codes_covered": 21679,
        "orbit_results": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        "PASS: new component/product exclusions and the Golay-factor obstruction replayed."
    )


if __name__ == "__main__":
    main()
