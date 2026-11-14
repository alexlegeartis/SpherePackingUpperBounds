"""Numerical necessary-contact test in catalog minimum-shell eigenspaces."""

from pathlib import Path
import json
import numpy as np
from scipy.optimize import linprog

P = Path(__file__).resolve().parent
O = P / "work" / "catalog"
out = []
for seed in json.loads((P / "catalog_screen.json").read_text()):
    if len(seed["moment_eigenvalues"]) == 1:
        continue
    name = seed["name"]
    data = np.load(O / f"{name}.npz")
    Q = data["Q"]
    V = np.loadtxt(O / f"{name}_vectors.txt", dtype=np.int64, ndmin=2)
    L = np.linalg.cholesky(Q)
    values, E = np.linalg.eigh(L.T @ (V.T @ V) @ L)
    groups = []
    for i, e in enumerate(values):
        if not groups or abs(e - values[groups[-1][0]]) > 1e-5:
            groups.append([i])
        else:
            groups[-1].append(i)
    parts = np.array([L @ E[:, g] @ E[:, g].T @ L.T for g in groups])
    A = np.einsum("vi,kij,vj->vk", V, parts, V)
    ranks = np.array(list(map(len, groups)))
    sol = linprog(ranks, A_ub=-A, b_ub=np.zeros(len(A)), bounds=(-1, 1), method="highs")
    row = {
        "name": name,
        "dimension": len(Q),
        "eigenspace_ranks": ranks.tolist(),
        "objective": float(sol.fun),
        "direction": sol.x.tolist(),
    }
    out.append(row)
    print(row, flush=True)
    np.savez_compressed(O / f"{name}_spectral.npz", parts=parts, ranks=ranks)
(P / "spectral_screen.json").write_text(json.dumps(out, indent=2) + "\n")
