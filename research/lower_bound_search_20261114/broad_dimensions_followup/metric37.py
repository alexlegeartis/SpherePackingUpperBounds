"""Full-metric discovery screen of the saved laminated 37D seed.

Root coordinates are in new36.npz's basis. This is a necessary-contact LP;
descent directions require further vectors and complete minimum verification.
"""

from pathlib import Path
import json, time, subprocess
import numpy as np
from scipy.optimize import linprog

P = Path(__file__).resolve().parent
work = P / "work"
work.mkdir(exist_ok=True)
seed = np.load(P / "laminated37.npz")
Q = seed["Q"]
n = len(Q)
if not (work / "new36_vectors.txt").exists():
    subprocess.run(
        [
            "g++",
            "-O3",
            "-std=c++17",
            str(P.parent / "experiments" / "enum_shell.cpp"),
            "-o",
            str(work / "enum_shell"),
        ],
        check=True,
    )
    np.savetxt(work / "new36_R.txt", np.linalg.cholesky(Q[:36, :36] * 20).T)
    subprocess.run(
        [
            str(work / "enum_shell"),
            "36",
            "160.0000001",
            "250000000",
            str(work / "new36_R.txt"),
            str(work / "new36_vectors.txt"),
            str(work / "new36_enum.log"),
        ],
        check=True,
    )
base = np.loadtxt(work / "new36_vectors.txt", dtype=np.int64, ndmin=2)
cross = np.loadtxt(P / "laminated37_cross.txt", dtype=np.int64, ndmin=2)
W = np.r_[
    np.c_[base, np.zeros(len(base), dtype=np.int64)],
    np.c_[cross, np.ones(len(cross), dtype=np.int64)],
]
assert np.max(abs(np.einsum("vi,ij,vj->v", W, Q, W) - 8)) < 1e-6
L = np.linalg.cholesky(Q)
X = W @ L
ii, jj = np.triu_indices(n)
c = (ii == jj).astype(float)
rng = np.random.default_rng(3737)
active = np.r_[rng.choice(len(base), 5000, replace=False), np.arange(len(base), len(W))]
rows = []
started = time.time()
for iteration in range(10):
    x = X[active]
    A = x[:, ii] * x[:, jj] * np.where(ii == jj, 1, 2)
    sol = linprog(
        c,
        A_ub=-A,
        b_ub=np.zeros(len(A)),
        bounds=(-1, 1),
        method="highs-ipm",
        options={
            "time_limit": 150,
            "run_crossover": "off",
            "ipm_optimality_tolerance": 1e-8,
        },
    )
    row = {
        "iteration": iteration,
        "constraints": len(active),
        "status": int(sol.status),
        "message": sol.message,
        "seconds": time.time() - started,
    }
    if sol.x is None:
        rows.append(row)
        print(row, flush=True)
        break
    D = np.zeros((n, n))
    D[ii, jj] = sol.x
    D[jj, ii] = sol.x
    residual = np.einsum("vi,ij,vj->v", X, D, X)
    bad = np.flatnonzero(residual < -1e-6)
    row.update(objective=float(sol.fun), violated=len(bad), worst=float(residual.min()))
    rows.append(row)
    print(row, flush=True)
    if sol.fun >= -1e-6 or not len(bad):
        np.savez_compressed(
            P / "metric37_direction.npz",
            D=D,
            L=L,
            Q=Q,
            active=active,
            weights=-sol.ineqlin.marginals,
        )
        break
    active = np.union1d(active, bad[np.argsort(residual[bad])[:4000]])
(P / "metric37_results.json").write_text(json.dumps(rows, indent=2) + "\n")
