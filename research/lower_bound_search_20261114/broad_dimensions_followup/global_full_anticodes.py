"""Extend the rank 2--4 minimum-root section searches to every beta < mu.

For a finite subset S of a torsion-free group, |S-S| >= 2|S|-1. The exact
dual-ball count therefore gives a useful cardinality upper bound before any
maximum-clique call. Unsolved node-capped graphs are explicitly retained.
"""

from pathlib import Path
import itertools, json, time, math, sys
from fractions import Fraction
import numpy as np
from flint import fmpz_mat
from dual_anticodes import maximum_clique

P = Path(__file__).resolve().parent
old = json.loads((P / "global_anticodes.json").read_text())
out = []
stem = "global_full_anticodes"
if len(sys.argv) > 1:
    mu = int(sys.argv[1])
    n = 19 if mu == 4 else 43
    stem = f"global_full_rank5_mu{mu}"
    r = json.loads(
        (
            P / ("rank5_results.json" if mu == 4 else "rank5_mu6_results.json")
        ).read_text()
    )
    old = [
        {
            "ambient_minimum": mu,
            "section_rank": 5,
            "dimension": n,
            "best": {
                "center_density_squared": r["best_squared"],
                "center_density": math.sqrt(float(Fraction(r["best_squared"]))),
            },
        }
    ]
    (P / f"section_Gram_types_mu{mu}_k5.json").write_text(
        json.dumps(np.load(P / f"section_Gram_types_mu{mu}_k5.npz")["types"].tolist())
        + "\n"
    )
for source in old:
    mu = source["ambient_minimum"]
    k = source["section_rank"]
    n = source["dimension"]
    types = json.loads((P / f"section_Gram_types_mu{mu}_k{k}.json").read_text())
    V = np.array(
        [v for v in itertools.product(range(1 - mu, mu), repeat=k) if any(v)],
        dtype=np.int64,
    )
    best = Fraction(source["best"]["center_density_squared"])
    bestrow = source["best"]
    graphs = 0
    unknown = []
    started = time.time()
    max_vertices = 0
    largestnodes = 0
    for typeid, raw in enumerate(types):
        K = np.array(raw, dtype=np.int64)
        F = fmpz_mat(raw)
        det = int(F.det())
        M = np.array((F.inv() * det).tolist(), dtype=np.int64)
        norm = np.einsum("vi,ij,vj->v", V, M, V)
        keep = (norm > 0) & (norm < mu * det)
        T = V[keep]
        norm = norm[keep]
        order = np.argsort(norm)
        T = T[order]
        norm = norm[order]
        minimum_checked = False
        for b in map(int, np.unique(norm)):
            j = int(np.searchsorted(norm, b, side="right"))
            upper_s = j // 2 + 1
            upper = Fraction(upper_s**2 * (mu * det - b) ** n, 4**n * det ** (n + 1))
            if upper <= best:
                continue
            X = T[:j]
            max_vertices = max(max_vertices, j)
            D = norm[:j, None] + norm[None, :j] - 2 * X @ M @ X.T
            A = D <= b
            np.fill_diagonal(A, False)
            adj = [sum(1 << int(v) for v in np.flatnonzero(row)) for row in A]
            try:
                clique, nodes = maximum_clique(adj, 1000000)
            except RuntimeError:
                unknown.append(
                    {
                        "type": typeid,
                        "diameter": str(Fraction(b, det)),
                        "vertices": j,
                        "cardinality_bound": upper_s,
                    }
                )
                continue
            graphs += 1
            largestnodes = max(largestnodes, nodes)
            s = len(clique) + 1
            score = Fraction(s * s * (mu * det - b) ** n, 4**n * det ** (n + 1))
            if score > best:
                if k == 5 and not minimum_checked:
                    bounds = [
                        math.isqrt((mu - 2) * int(M[i, i]) // det) for i in range(k)
                    ]
                    small = np.array(
                        list(itertools.product(*(range(-v, v + 1) for v in bounds))),
                        dtype=np.int64,
                    )
                    vals = np.einsum("vi,ij,vj->v", small, K, small)
                    minimum_checked = True
                    if np.any((vals > 0) & (vals < mu)):
                        break
                best = score
                bestrow = {
                    "K": raw,
                    "determinant": det,
                    "diameter": str(Fraction(b, det)),
                    "cluster_size": s,
                    "cluster": [[0] * k] + X[clique].tolist(),
                    "center_density_squared": str(score),
                    "center_density": math.sqrt(float(score)),
                }
                print("IMPROVEMENT", n, bestrow, flush=True)
        if typeid % 100 == 0:
            print(
                "progress",
                n,
                typeid,
                "graphs",
                graphs,
                "unknown",
                len(unknown),
                "max_vertices",
                max_vertices,
                "seconds",
                time.time() - started,
                flush=True,
            )
    row = {
        "dimension": n,
        "ambient_minimum": mu,
        "section_rank": k,
        "diameter_range": f"0 <= beta < {mu}",
        "Gram_types": len(types),
        "graphs_solved": graphs,
        "maximum_graph_vertices": max_vertices,
        "maximum_clique_search_nodes": largestnodes,
        "unknown_graphs": unknown,
        "best": bestrow,
        "seconds": time.time() - started,
    }
    out.append(row)
    print(
        "DONE",
        n,
        "best",
        bestrow["center_density"],
        "unknown",
        len(unknown),
        flush=True,
    )
    (P / f"{stem}.json").write_text(json.dumps(out, indent=2) + "\n")
