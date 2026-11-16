"""Deterministic overlattice discovery screen. Requires fpylll; not a verifier."""

from pathlib import Path
import numpy as np, json, time
from flint import fmpz_mat, fmpq_mat
from fpylll import IntegerMatrix, GSO, LLL, BKZ, FPLLL
import argparse

P = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser(
    description="Finite discovery only; all minimum claims require separate verification."
)
ap.add_argument("--output", type=Path, required=True)
args = ap.parse_args()
if args.output.exists():
    raise SystemExit("Output already exists; choose a new path.")
FPLLL.set_precision(160)
rng = np.random.default_rng(202611165456)
out = []
start = time.time()
for n, prime, index, trials in [(54, 3, 9, 24), (56, 29, 29, 24)]:
    d = np.load(
        P
        / (
            "cm58_20261116/frontier/local54_55.npz"
            if n == 54
            else "cm58_20261116/certificate.npz"
        )
    )
    G = fmpz_mat(d["G54" if n == 54 else "G"].tolist())
    dual = prime * G.inv()
    assert all(x.denominator == 1 for row in dual.tolist() for x in row)
    D = np.array(
        [[int(x) % prime for x in row] for row in dual.tolist()], dtype=np.int64
    )
    best = 0
    for trial in range(trials):
        r = rng.integers(0, prime, n, dtype=np.int64) @ D % prime
        generators = (prime * np.eye(n, dtype=np.int64)).tolist() + [r.tolist()]
        if n == 54:
            generators.append(d["q3"][0].tolist())
        B = fmpz_mat(fmpz_mat(generators).hnf().tolist()[:n])
        if abs(B.det()) != prime**n // index:
            continue
        full = B * G * B.transpose()
        assert all(int(x) % prime == 0 for row in full.tolist() for x in row)
        Q = fmpz_mat([[int(x) // prime for x in row] for row in full.tolist()])
        A = IntegerMatrix.from_matrix([[int(x) for x in row] for row in Q.tolist()])
        U = IntegerMatrix.identity(n)
        M = GSO.Mat(A, U=U, gram=True, float_type="mpfr", update=True)
        lll = LLL.Reduction(M)
        lll()
        BKZ.Reduction(M, lll, BKZ.Param(block_size=28, max_loops=2))()
        T = fmpz_mat([[int(U[i, j]) for j in range(n)] for i in range(n)])
        R = T * Q * T.transpose()
        at = min(range(n), key=lambda i: R[i, i])
        norm = int(R[at, at])
        w = list(map(int, (T * B).tolist()[at]))
        assert (
            int((fmpz_mat([w]) * G * fmpz_mat([w]).transpose())[0, 0]) == prime * norm
        )
        best = max(best, norm / prime)
        ratio = index * (norm / (10 * prime)) ** (n / 2)
        out.append(
            {
                "dimension": n,
                "trial": trial,
                "prime": prime,
                "index": index,
                "added_numerator": r.tolist(),
                "witness_numerator": w,
                "witness_squared_norm": f"{norm}/{prime}",
                "density_ratio_upper_float": ratio,
                "rejected": ratio <= 1,
            }
        )
        args.output.write_text(json.dumps(out, indent=2) + "\n")
        print(
            n,
            trial,
            "norm upper",
            norm / prime,
            "best",
            best,
            "density ratio upper",
            ratio,
            "seconds",
            time.time() - start,
            flush=True,
        )
