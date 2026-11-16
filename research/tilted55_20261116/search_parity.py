"""Deterministic finite discovery screen; not a sphere-packing verifier.

Requires NumPy, python-flint and fpylll. It generates supplied orbits anew,
then tries 20 local-search restarts per dimension. Every positive minimum
claim would require a separate complete certificate.
"""

from pathlib import Path
import numpy as np, json, time
from flint import nmod_mat, fmpz_mat
from fpylll import IntegerMatrix, GSO, LLL, BKZ, FPLLL
import argparse, sys

P = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(P / "cm58_20261116"))
import verify_frontier as frontier

ap = argparse.ArgumentParser(
    description="Capped discovery only; no positive packing certificates."
)
ap.add_argument("--output", type=Path, required=True)
args = ap.parse_args()
if args.output.exists():
    raise SystemExit("Output already exists; choose a new path.")
FPLLL.set_precision(160)
rng = np.random.default_rng(202611165480)
d54 = np.load(P / "cm58_20261116/frontier/local54_55.npz", allow_pickle=False)
d56 = np.load(P / "cm58_20261116/certificate.npz", allow_pickle=False)
g56 = fmpz_mat(d56["G"].tolist())
v56, _ = frontier.orbit(
    g56,
    d56["point_coordinates"][0],
    np.load(P / "cm58_20261116/frontier/isometries56.npz")["maps"],
    10,
)
g80, seed80 = frontier.parent80()
v80, _ = frontier.orbit(
    g80, seed80, np.load(P / "cm58_20261116/frontier/isometries80.npz")["maps"], 14
)
sources = [
    (54, d54["G54"], d54["roots54"], 10),
    (56, d56["G"], v56, 10),
    (80, np.array(g80.tolist(), dtype=np.int64), v80, 14),
]
report = []
for n, G, V, minimum in sources:
    bits = V % 2
    best = len(V)
    found = 0
    start = time.time()
    probes = []
    for restart in range(20):
        H = rng.integers(0, 2, (7, n), dtype=np.int64)
        synd = (bits @ H.T % 2) @ (1 << np.arange(7))
        stall = 0
        for step in range(2000):
            loss = int(np.count_nonzero(synd == 0))
            best = min(best, loss)
            if loss == 0:
                break
            zero = bits[synd == 0].sum(axis=0)
            gain = np.array(
                [zero - bits[synd == (1 << j)].sum(axis=0) for j in range(7)]
            )
            maximum = int(gain.max())
            choices = np.argwhere(gain == maximum)
            j, k = map(int, choices[int(rng.integers(len(choices)))])
            if maximum <= 0:
                stall += 1
                if stall % 30 == 0:
                    j = int(rng.integers(7))
                    k = int(rng.integers(n))
            else:
                stall = 0
            H[j, k] ^= 1
            synd[bits[:, k].astype(bool)] ^= 1 << j
        if loss == 0 and nmod_mat(H.tolist(), 2).rank() == 7:
            found += 1
            N, k = nmod_mat(H.tolist(), 2).nullspace()
            basis = [[int(N[i, j]) for i in range(n)] for j in range(k)] + (
                2 * np.eye(n, dtype=int)
            ).tolist()
            B = fmpz_mat(fmpz_mat(basis).hnf().tolist()[:n])
            assert abs(B.det()) == 128
            Q = B * fmpz_mat(G.tolist()) * B.transpose()
            q = IntegerMatrix.from_matrix([[int(x) for x in row] for row in Q.tolist()])
            U = IntegerMatrix.identity(n)
            M = GSO.Mat(q, U=U, gram=True, float_type="mpfr", update=True)
            lll = LLL.Reduction(M)
            lll()
            BKZ.Reduction(M, lll, BKZ.Param(block_size=26, max_loops=2))()
            T = fmpz_mat([[int(U[i, j]) for j in range(n)] for i in range(n)])
            R = T * Q * T.transpose()
            at = min(range(n), key=lambda i: R[i, i])
            w = list(map(int, (T * B).tolist()[at]))
            norm = int(R[at, at])
            assert (fmpz_mat([w]) * fmpz_mat(G.tolist()) * fmpz_mat([w]).transpose())[
                0, 0
            ] == norm
            assert np.all(np.array(w) @ H.T % 2 == 0)
            probes.append(
                {
                    "checks": H.tolist(),
                    "known_rays_deleted": len(V),
                    "found_norm_upper": norm,
                    "counterexample": list(map(int, w)),
                    "rejected": norm == minimum,
                }
            )
            print(
                "PROBE",
                n,
                restart,
                "known",
                len(V),
                "norm",
                norm,
                "seconds",
                time.time() - start,
                flush=True,
            )
            if found >= 5:
                break
        if time.time() - start > 90:
            break
    report.append(
        {
            "dimension": n,
            "index": 128,
            "known_rays": len(V),
            "best_known_shell_survivors": best,
            "candidates": probes,
            "seconds": time.time() - start,
            "record_claimed": False,
        }
    )
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        "DONE",
        n,
        "best",
        best,
        "found",
        found,
        "seconds",
        time.time() - start,
        flush=True,
    )
