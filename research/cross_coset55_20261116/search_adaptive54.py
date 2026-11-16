"""Adaptive index-nine overlattice discovery; requires fpylll."""

from pathlib import Path
import numpy as np, json, time
from flint import fmpz_mat
from fpylll import IntegerMatrix, GSO, LLL, BKZ, FPLLL
import argparse

P = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser(
    description="Adaptive finite discovery only; no positive certificate."
)
ap.add_argument("--output", type=Path, required=True)
ap.add_argument("--trials", type=int, default=300)
args = ap.parse_args()
if args.output.exists():
    raise SystemExit("Output already exists; choose a new path.")
FPLLL.set_precision(160)
rng = np.random.default_rng(545600)
d = np.load(P / "cm58_20261116/frontier/local54_55.npz")
G = fmpz_mat(d["G54"].tolist())
D = np.array(
    [[int(x) % 3 for x in row] for row in (3 * G.inv()).tolist()], dtype=np.int64
)
initial = max(
    [
        a
        for a in json.loads(
            (P / "cross_coset55_20261116/overlattice_screen.json").read_text()
        )["rows"]
        if a["dimension"] == 54
    ],
    key=lambda a: float(__import__("fractions").Fraction(a["witness_squared_norm"])),
)
current = np.array(initial["added_numerator"])
best = 0
score = ()
records = []
start = time.time()
for step in range(args.trials):
    r = (current + int(rng.integers(1, 3)) * D[int(rng.integers(54))]) % 3
    if step % 50 == 49:
        r = rng.integers(0, 3, 54) @ D % 3
        score = ()
    generators = (3 * np.eye(54, dtype=np.int64)).tolist() + [
        r.tolist(),
        d["q3"][0].tolist(),
    ]
    B = fmpz_mat(fmpz_mat(generators).hnf().tolist()[:54])
    if abs(B.det()) != 3**52:
        continue
    Q = B * G * B.transpose()
    Q = fmpz_mat([[int(x) // 3 for x in row] for row in Q.tolist()])
    A = IntegerMatrix.from_matrix([[int(x) for x in row] for row in Q.tolist()])
    U = IntegerMatrix.identity(54)
    M = GSO.Mat(A, U=U, gram=True, float_type="mpfr", update=True)
    lll = LLL.Reduction(M)
    lll()
    BKZ.Reduction(M, lll, BKZ.Param(block_size=32, max_loops=2))()
    T = fmpz_mat([[int(U[i, j]) for j in range(54)] for i in range(54)])
    R = T * Q * T.transpose()
    norms = sorted(int(R[i, i]) for i in range(54))
    newscore = tuple(norms[:12])
    at = min(range(54), key=lambda i: R[i, i])
    w = list(map(int, (T * B).tolist()[at]))
    norm = norms[0]
    records.append(
        {
            "step": step,
            "added_numerator": r.tolist(),
            "norm_numerator": norm,
            "norm_denominator": 3,
            "witness_numerator": w,
        }
    )
    if newscore >= score:
        current = r
        score = newscore
    if norm > best:
        best = norm
        print("BEST", step, norm, "/3", newscore, flush=True)
    args.output.write_text(json.dumps(records, indent=2) + "\n")
    if step % 25 == 0:
        print(
            "STEP",
            step,
            "best",
            best,
            "current",
            norm,
            "seconds",
            time.time() - start,
            flush=True,
        )
    if norm >= 28:
        print("TARGET_CANDIDATE_REQUIRES_CERTIFICATE", flush=True)
        break
