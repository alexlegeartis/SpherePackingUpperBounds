"""Replay the small exact searches, contact certificate and rejection witnesses.

--full also repeats the longer 19D and 43D rank-five searches. All generated
files and executables stay in a temporary directory.
"""

from pathlib import Path
import itertools, json, math, shutil, subprocess, sys, tempfile
from fractions import Fraction
import numpy as np
from dual_anticodes import maximum_clique

P = Path(__file__).resolve().parent


def brute_clique(adj):
    best = 0
    for mask in range(1 << len(adj)):
        if mask.bit_count() <= best:
            continue
        if all(
            (mask & ~(1 << i) & ~adj[i]) == 0 for i in range(len(adj)) if mask >> i & 1
        ):
            best = mask.bit_count()
    return best


for n in range(6):
    edges = list(itertools.combinations(range(n), 2))
    for code in range(1 << len(edges)):
        adj = [0] * n
        for j, (a, b) in enumerate(edges):
            if code >> j & 1:
                adj[a] |= 1 << b
                adj[b] |= 1 << a
        clique, _ = maximum_clique(adj)
        assert len(clique) == brute_clique(adj)
rng = np.random.default_rng(20261114)
for _ in range(50):
    adj = [0] * 11
    for a, b in itertools.combinations(range(11), 2):
        if rng.integers(2):
            adj[a] |= 1 << b
            adj[b] |= 1 << a
    assert len(maximum_clique(adj)[0]) == brute_clique(adj)
print("PASS: clique search versus independent exhaustive subset checks", flush=True)


def orbit_count(base, k=5):
    edges = list(itertools.combinations(range(k), 2))
    lookup = {e: i for i, e in enumerate(edges)}
    fixed = 0
    group = 0
    for perm in itertools.permutations(range(k)):
        for signs0 in itertools.product((-1, 1), repeat=k - 1):
            signs = (1,) + signs0
            image = [lookup[tuple(sorted((perm[i], perm[j])))] for i, j in edges]
            sgn = [signs[i] * signs[j] for i, j in edges]
            seen = set()
            positive = 0
            for start in range(len(edges)):
                if start in seen:
                    continue
                at = start
                product = 1
                while at not in seen:
                    seen.add(at)
                    product *= sgn[at]
                    at = image[at]
                positive += product == 1
            fixed += base**positive
            group += 1
    assert fixed % group == 0
    return fixed // group


assert orbit_count(5) == 6595 and orbit_count(7) == 160419
print("PASS: independent Burnside counts for rank-five Gram classes", flush=True)

repo = P.parents[2]
data = repo / "constructions" / "dim36" / "data"
U = np.loadtxt(data / "U.txt", dtype=np.int64)
original = np.loadtxt(data / "Q_original.txt", dtype=np.int64)
Q = U @ original @ U.T
for q, cut, count in [(3, 1120, 47444), (4, 2000, 44831)]:
    d = np.load(P / f"deep_hole_witnesses_q{q}.npz")
    Y = d["witnesses"].astype(np.int64)
    residues = d["residues"].astype(np.int64)
    assert len(Y) == count and len(np.unique(d["ids"])) == count
    assert np.all((Y - residues) % q == 0)
    assert 36**2 * int(abs(Q).max()) * int(abs(Y).max()) ** 2 < 2**63
    norm = np.einsum("vi,ij,vj->v", Y, Q, Y)
    assert np.array_equal(norm, d["norms"])
    assert np.all(norm <= cut)
    print(f"PASS: {count} exact rejecting witnesses for denominator {q}", flush=True)

known = {
    19: Fraction(13**19, 3**18 * 5**21),
    20: Fraction(3**40, 2**20 * 5**21),
    21: Fraction(43**21, 2**82 * 3**23),
    22: Fraction(11**22, 2**46 * 3**21),
    43: Fraction(137**43, 2**223 * 3**42),
    44: Fraction(157**44, 2**44 * 5**43 * 11**46),
    45: Fraction(23**45, 2**183),
    46: Fraction(13**46, 3**93),
}
for N, mu in [(24, 4), (48, 6)]:
    one = Fraction(mu, 4) ** (N - 1) / mu
    two = Fraction(4 * mu, N) * Fraction(mu * (N - 1), 4 * N) ** (N - 1)
    assert two > one and N % mu == 0
print("PASS: exact rank-one optimum comparisons for 23D and 47D", flush=True)

with tempfile.TemporaryDirectory(prefix="packing-broad-check-") as directory:
    work = Path(directory)
    for pattern in ("*.py", "*.npz", "*.json", "*.cpp"):
        for file in P.glob(pattern):
            shutil.copy2(file, work / file.name)

    def run(*args):
        result = subprocess.run(
            args, cwd=work, check=True, capture_output=True, text=True
        )
        return result.stdout

    for mu, expected in [
        (4, "classes 6595 relaxed_types 2006"),
        (6, "classes 160419 relaxed_types 51647"),
    ]:
        if mu == 4:
            run(
                "g++",
                "-O3",
                "-std=c++17",
                "gram_types5.cpp",
                "-o",
                str(work / "gram_types5"),
            )
        generated = work / f"types{mu}.txt"
        summary = run(str(work / "gram_types5"), str(generated), str(mu)).strip()
        assert summary == expected
        assert np.array_equal(
            np.loadtxt(generated, dtype=np.int64, ndmin=2).reshape(-1, 5, 5),
            np.load(work / f"section_Gram_types_mu{mu}_k5.npz")["types"],
        )
    print("PASS: exact rank-five class generation matches archived inputs", flush=True)
    run(sys.executable, "global_anticodes.py")
    run(sys.executable, "global_full_anticodes.py")
    result = json.loads((work / "global_full_anticodes.json").read_text())
    assert len(result) == 6 and all(not row["unknown_graphs"] for row in result)
    bounded = json.loads((work / "global_anticodes.json").read_text())
    for a, b in zip(result, bounded):
        assert (
            Fraction(a["best"]["center_density_squared"])
            == Fraction(b["best"]["center_density_squared"])
            == known[a["dimension"]]
        )
    print("PASS: all-diameter rank 2--4 searches in 20--22D and 44--46D", flush=True)
    run(sys.executable, "rank2_long_sections.py")
    longer = json.loads((work / "rank2_long_sections.json").read_text())
    assert all(
        not r["unknown_graphs"] and Fraction(r["best_squared"]) == known[r["dimension"]]
        for r in longer
    )
    print(
        "PASS: all rank-two even reduced section Grams through basis norm 64",
        flush=True,
    )
    print(run(sys.executable, "certify_periodic44.py"), end="", flush=True)
    if "--full" in sys.argv:
        run(sys.executable, "rank5_anticodes.py", "4")
        run(sys.executable, "global_full_anticodes.py", "4")
        r = json.loads((work / "rank5_results.json").read_text())
        assert (
            Fraction(r["best_squared"]) == Fraction(r["record_squared"])
            and not r["unknown_graphs"]
        )
        full = json.loads((work / "global_full_rank5_mu4.json").read_text())[0]
        assert (
            Fraction(full["best"]["center_density_squared"])
            == Fraction(r["record_squared"])
            and not full["unknown_graphs"]
        )
        print("PASS: full 19D rank-five search", flush=True)
        run(sys.executable, "rank5_anticodes.py", "6")
        r = json.loads((work / "rank5_mu6_results.json").read_text())
        assert (
            Fraction(r["best_squared"]) == Fraction(r["record_squared"])
            and not r["unknown_graphs"]
        )
        print("PASS: bounded-diameter 43D rank-five search", flush=True)
print("PASS: broad-dimension follow-up checks; no new SOTA asserted", flush=True)
