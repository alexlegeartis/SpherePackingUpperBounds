"""All even reduced rank-two section Grams with basis norms at most 64."""

from pathlib import Path
import math, json, time
from fractions import Fraction
import numpy as np
from dual_anticodes import maximum_clique

P = Path(__file__).resolve().parent
out = []
for mu, N, record in [
    (4, 24, Fraction(3, 4) * Fraction(11, 12) ** 22),
    (6, 48, Fraction(1, 3) * Fraction(13, 9) ** 46),
]:
    n = N - 2
    best = record
    winner = None
    graphs = 0
    types = 0
    unknown = []
    started = time.time()
    for a in range(mu, 65, 2):
        for c in range(a, 65, 2):
            V = np.array(
                [
                    (x, y)
                    for x in range(-math.isqrt(mu * a - 1), math.isqrt(mu * a - 1) + 1)
                    for y in range(-math.isqrt(mu * c - 1), math.isqrt(mu * c - 1) + 1)
                    if x or y
                ],
                dtype=np.int64,
            )
            for b in range(a // 2 + 1):
                types += 1
                det = a * c - b * b
                M = np.array([[c, -b], [-b, a]])
                norm = np.einsum("vi,ij,vj->v", V, M, V)
                keep = norm < mu * det
                T = V[keep]
                norm = norm[keep]
                order = np.argsort(norm)
                T = T[order]
                norm = norm[order]
                assert Fraction(mu, 4) ** n / det <= record
                for t in map(int, np.unique(norm)):
                    j = int(np.searchsorted(norm, t, side="right"))
                    upper = Fraction(
                        (j // 2 + 1) ** 2 * (mu * det - t) ** n, 4**n * det ** (n + 1)
                    )
                    if upper <= best:
                        continue
                    X = T[:j]
                    D = norm[:j, None] + norm[None, :j] - 2 * X @ M @ X.T
                    A = D <= t
                    np.fill_diagonal(A, False)
                    adj = [sum(1 << int(v) for v in np.flatnonzero(row)) for row in A]
                    try:
                        clique, nodes = maximum_clique(adj, 1000000)
                    except RuntimeError:
                        unknown.append([a, b, c, t])
                        continue
                    graphs += 1
                    s = len(clique) + 1
                    score = Fraction(s * s * (mu * det - t) ** n, 4**n * det ** (n + 1))
                    if score > best:
                        best = score
                        winner = {
                            "K": [[a, b], [b, c]],
                            "cluster": [[0, 0]] + X[clique].tolist(),
                            "diameter": str(Fraction(t, det)),
                        }
                        print(
                            "IMPROVEMENT",
                            n,
                            winner,
                            math.sqrt(float(score)),
                            flush=True,
                        )
        print(
            "progress",
            n,
            a,
            "types",
            types,
            "graphs",
            graphs,
            "seconds",
            time.time() - started,
            flush=True,
        )
    out.append(
        {
            "dimension": n,
            "max_basis_norm": 64,
            "Gram_types": types,
            "graphs_solved": graphs,
            "unknown_graphs": unknown,
            "best_squared": str(best),
            "record_squared": str(record),
            "winner": winner,
            "seconds": time.time() - started,
        }
    )
    (P / "rank2_long_sections.json").write_text(json.dumps(out, indent=2) + "\n")
