"""Heuristic coordinate-permutation search across published ternary code types."""

from pathlib import Path
import argparse
import json
import subprocess
import numpy as np
from finite import (
    parent_supports,
    shortened_supports,
    write_instance,
    ternary_code,
    audit_selforthogonal,
)

P = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dimension", type=int, choices=(30, 32, 34, 36, 38, 40), required=True
    )
    parser.add_argument("--steps", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261114)
    parser.add_argument("--output-prefix", type=Path)
    args = parser.parse_args()
    n = args.dimension // 2
    codes = np.load(P / "inputs/catalogs.npz")[str(n)]
    assert len(codes) == {15: 12, 16: 7, 17: 23, 18: 160, 19: 56, 20: 24}[n]
    for G in codes:
        assert not (G @ G.T % 3).any()
        ternary_code(G)  # Checks rank by complete finite enumeration.
    prefix = (args.output_prefix or P / f"work/catalog{n}").resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    supports, _, _ = parent_supports()
    current = shortened_supports(supports, list(range(n, 20)))
    instance = Path(str(prefix) + "_instance.txt")
    catalog = Path(str(prefix) + "_catalog.txt")
    binary = Path(str(prefix) + "_search")
    write_instance(instance, n, codes.shape[1], current)
    catalog.write_text(
        f"{len(codes)} {n}\n"
        + "\n".join(" ".join(map(str, row)) for G in codes for row in G)
        + "\n"
    )
    subprocess.run(
        [
            "g++",
            "-std=c++17",
            "-O3",
            str(P / "src/permutation_search.cpp"),
            "-o",
            str(binary),
        ],
        check=True,
    )
    with Path(str(prefix) + ".log").open("w") as log:
        subprocess.run(
            [
                str(binary),
                str(instance),
                str(catalog),
                str(prefix),
                str(args.steps),
                str(args.seed),
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    best = Path(str(prefix) + ".best.txt")
    lines = best.read_text().splitlines()
    header = list(map(int, lines[0].split()))
    G = np.array([list(map(int, line.split())) for line in lines[1:]], dtype=np.int16)
    result = audit_selforthogonal(G, current)
    assert result["search_energy"] == header[2]
    result.update(
        catalog_types=len(codes),
        steps_per_type=args.steps,
        seed=args.seed,
        scope="Heuristic placements; not exhaustive over coordinate permutations.",
    )
    Path(str(prefix) + "_audit.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
