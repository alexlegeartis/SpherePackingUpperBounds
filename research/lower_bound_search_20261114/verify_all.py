"""Replay the exact spectral calculation and exhaustive fiber obstruction."""
from pathlib import Path
import itertools
import subprocess
import sys
import tempfile

P = Path(__file__).resolve().parent
assert sys.flags.optimize == 0, "Do not disable assertions"


def run(command):
    result = subprocess.run([str(x) for x in command], check=True,
                            text=True, capture_output=True)
    print(result.stdout, end="", flush=True)
    return result.stdout


def audit_sat(text, n, supports):
    lines = text.splitlines()
    assert lines[0].startswith("SAT nodes ")
    columns = [tuple(map(int, x.split())) for x in lines[1:]]
    assert len(columns) == n and all(len(x) == 3 and any(x) for x in columns)
    for support in supports:
        indices = [i for i in range(n) if support >> i & 1]
        # Independent direct enumeration, without bit permutations or search state.
        for signs in itertools.product((-1, 1), repeat=len(indices)):
            value = tuple(sum(sign * columns[i][j] for i, sign in zip(indices, signs)) % 3
                          for j in range(3))
            assert value != (0, 0, 0)


with tempfile.TemporaryDirectory(prefix="lower-bound-replay-") as directory:
    W = Path(directory)
    run([sys.executable, P / "verify_codes.py", "--output-dir", W])
    # A positive instance must be accepted and its answer checked independently.
    (W / "positive.txt").write_text("6 1\n63\n")
    for name, flags in [("fast", ["-O3"]),
                        ("coordinate_check", ["-O2", "-DINDEPENDENT_STEPS"])]:
        binary = W / name
        run(["g++", "-std=c++17", *flags, P / "fiber_csp.cpp", "-o", binary])
        audit_sat(run([binary, W / "positive.txt", "1000000"]), 6, [63])
        for n in (17, 18):
            text = run([binary, W / f"supports{n}.txt", "1000000"])
            assert text.startswith("UNSAT nodes 23687 "), text
        if name == "fast":
            for n in (17, 18):
                text = run([binary, W / f"supports{n}.txt", "20000000", "full"])
                assert text.startswith("UNSAT nodes 1683150 "), text
    print("PASS: both finite obstructions replayed, including searches without symmetry reduction.")
