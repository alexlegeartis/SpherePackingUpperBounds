"""Exact audit of a saved self-orthogonal minimum-twelve fiber candidate."""

from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
from finite import parent_supports, shortened_supports, audit_selforthogonal

P = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--removed", type=int, nargs="*")
    parser.add_argument("--output", type=Path, default=P / "work/audit.json")
    args = parser.parse_args()
    lines = args.candidate.read_text().splitlines()
    header = list(map(int, lines[0].split()))
    n, t = header[:2]
    assert n in range(15, 21)
    G = np.array([list(map(int, line.split())) for line in lines[1:]], dtype=np.int16)
    assert G.shape == (t, n)
    selected = args.removed if args.removed is not None else list(range(n, 20))
    assert len(selected) == 20 - n and len(set(selected)) == len(selected)
    supports, _, _ = parent_supports()
    result = audit_selforthogonal(G, shortened_supports(supports, selected))
    penalty = header[4] if len(header) >= 5 else 10
    expected = (
        penalty * result["norm6_bad_trit_words"] + result["norm9_bad_trit_words"]
    ) // 2
    assert expected == header[2]
    result.update(
        removed_coordinates=selected,
        header_norm6_penalty=penalty,
        ternary_generator_sha256=hashlib.sha256(
            G.astype(np.uint8).tobytes()
        ).hexdigest(),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
