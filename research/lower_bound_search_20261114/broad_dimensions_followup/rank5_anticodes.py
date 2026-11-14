"""Exact rank-five relaxed Gram/anticode search for 19D; no embedding assumed."""

from pathlib import Path
import json, time, math, sys
from fractions import Fraction
import numpy as np
from flint import fmpz_mat
from dual_anticodes import maximum_clique

P = Path(__file__).resolve().parent
mu = int(sys.argv[1]) if len(sys.argv) > 1 else 4
n = 19 if mu == 4 else 43
limit = Fraction(99, 100) if mu == 4 else Fraction(13, 20)
stem = "rank5" if mu == 4 else "rank5_mu6"
types = np.load(P / f"section_Gram_types_mu{mu}_k5.npz")["types"]
V = np.array(
    [v for v in __import__("itertools").product((-1, 0, 1), repeat=5) if any(v)],
    dtype=np.int64,
)
best = (
    Fraction(36, 300) * Fraction(13, 15) ** 19
    if mu == 4
    else Fraction(36, 3072) * Fraction(137, 96) ** 43
)
record = best
improvements = []
unknown = []
invalid = []
started = time.time()
graphs = 0
(P / f"{stem}_improvements.json").write_text("[]\n")


def short_witness(K, M, det):
    bounds = [math.isqrt((mu - 2) * int(M[i, i]) // det) for i in range(5)]
    T = np.array(
        list(__import__("itertools").product(*(range(-b, b + 1) for b in bounds))),
        dtype=np.int64,
    )
    norm = np.einsum("vi,ij,vj->v", T, K, T)
    bad = np.flatnonzero((norm > 0) & (norm < mu))
    return None if not len(bad) else (T[bad[0]].tolist(), int(norm[bad[0]]))


for typeid, K in enumerate(types):
    F = fmpz_mat(K.tolist())
    det = int(F.det())
    M = np.array((F.inv() * det).tolist(), dtype=np.int64)
    norm = np.einsum("vi,ij,vj->v", V, M, V)
    keep = (norm > 0) & (limit.denominator * norm <= limit.numerator * det)
    T = V[keep]
    norm = norm[keep]
    order = np.argsort(norm)
    T = T[order]
    norm = norm[order]
    minimum_checked = False
    for b in [0] + list(map(int, np.unique(norm))):
        j = int(np.searchsorted(norm, b, side="right"))
        upper = Fraction((j // 2 + 1) ** 2 * (mu * det - b) ** n, 4**n * det ** (n + 1))
        if upper <= best:
            continue
        X = T[:j]
        D = norm[:j, None] + norm[None, :j] - 2 * X @ M @ X.T
        A = D <= b
        np.fill_diagonal(A, False)
        adj = [sum(1 << int(v) for v in np.flatnonzero(row)) for row in A]
        try:
            clique, nodes = maximum_clique(adj, 1000000)
        except RuntimeError:
            unknown.append([typeid, str(Fraction(b, det))])
            continue
        graphs += 1
        s = len(clique) + 1
        score = Fraction(s * s * (mu * det - b) ** n, 4**n * det ** (n + 1))
        if score > best:
            if not minimum_checked:
                witness = short_witness(K, M, det)
                minimum_checked = True
                if witness is not None:
                    invalid.append(
                        {"type": typeid, "vector": witness[0], "norm": witness[1]}
                    )
                    break
            best = score
            row = {
                "type": typeid,
                "K": K.tolist(),
                "determinant": det,
                "diameter": str(Fraction(b, det)),
                "cluster_size": s,
                "cluster": [[0] * 5] + X[clique].tolist(),
                "score_squared": str(score),
                "density": math.sqrt(float(score)),
                "record_ratio": math.sqrt(float(score / record)),
            }
            improvements.append(row)
            print("IMPROVEMENT", json.dumps(row), flush=True)
            (P / f"{stem}_improvements.json").write_text(
                json.dumps(improvements, indent=2) + "\n"
            )
    if typeid % 100 == 0:
        print(
            "progress",
            typeid,
            "graphs",
            graphs,
            "unknown",
            len(unknown),
            "seconds",
            time.time() - started,
            flush=True,
        )
(P / f"{stem}_results.json").write_text(
    json.dumps(
        {
            "dimension": n,
            "Gram_types": len(types),
            "diameter_limit": str(limit),
            "graphs_solved": graphs,
            "unknown_graphs": unknown,
            "invalid_Gram_candidates": invalid,
            "best_squared": str(best),
            "record_squared": str(record),
            "improvements": improvements,
            "seconds": time.time() - started,
        },
        indent=2,
    )
    + "\n"
)
