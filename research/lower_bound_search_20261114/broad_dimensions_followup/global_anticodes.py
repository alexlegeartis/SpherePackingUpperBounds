"""Exact finite search of minimum-root section Grams and dual anticodes.

No embedding assumption is used. Matrices that pass the finite short-vector
filter form a SUPERSET of possible section Grams. Thus a failure to beat the
record excludes every realizable Gram in this explicitly bounded family.
"""

from pathlib import Path
import itertools, json, math, time
from fractions import Fraction
import numpy as np
from flint import fmpz_mat
from dual_anticodes import maximum_clique

P = Path(__file__).resolve().parent
out = []
for mu, N, k, limit in [
    (4, 24, 2, Fraction(99, 100)),
    (4, 24, 3, Fraction(99, 100)),
    (4, 24, 4, Fraction(99, 100)),
    (6, 48, 2, Fraction(13, 20)),
    (6, 48, 3, Fraction(13, 20)),
    (6, 48, 4, Fraction(13, 20)),
]:
    n = N - k
    ij = np.triu_indices(k, 1)
    perms = list(itertools.permutations(range(k)))
    signs = [np.array((1,) + s) for s in itertools.product((-1, 1), repeat=k - 1)]
    seen = set()
    types = []
    V = np.array(
        [v for v in itertools.product((-1, 0, 1), repeat=k) if any(v)], dtype=np.int64
    )
    started = time.time()
    all_grams = 0
    assert limit * mu < 4
    for key in itertools.product(range(-mu // 2, mu // 2 + 1), repeat=k * (k - 1) // 2):
        if key in seen:
            continue
        K = np.eye(k, dtype=np.int64) * mu
        K[ij] = key
        K[(ij[1], ij[0])] = key
        for perm in perms:
            B = K[np.ix_(perm, perm)]
            for s in signs:
                seen.add(tuple(map(int, (s[:, None] * B * s[None, :])[ij])))
        all_grams += 1
        if np.any(np.einsum("vi,ij,vj->v", V, K, V) < mu):
            continue
        if any(fmpz_mat(K[:j, :j].tolist()).det() <= 0 for j in range(1, k + 1)):
            continue
        types.append(K)
    (P / f"section_Gram_types_mu{mu}_k{k}.json").write_text(
        json.dumps([K.tolist() for K in types]) + "\n"
    )
    best = Fraction(0)
    bestrow = None
    graphs = 0
    for typeid, K in enumerate(types):
        F = fmpz_mat(K.tolist())
        det = int(F.det())
        M = np.array((F.inv() * det).tolist(), dtype=np.int64)
        assert np.array_equal(K @ M, det * np.eye(k, dtype=np.int64))
        norm = np.einsum("vi,ij,vj->v", V, M, V)
        keep = limit.denominator * norm <= limit.numerator * det
        T = V[keep]
        norm = norm[keep]
        D = norm[:, None] + norm[None, :] - 2 * T @ M @ T.T
        for b in [0] + sorted(set(map(int, norm))):
            ix = np.flatnonzero(norm <= b)
            upper = Fraction(
                (len(ix) + 1) ** 2 * (mu * det - b) ** n, 4**n * det ** (n + 1)
            )
            if upper < best:
                continue
            if not len(ix):
                clique = []
                nodes = 0
            else:
                A = D[np.ix_(ix, ix)] <= b
                np.fill_diagonal(A, False)
                adj = [sum(1 << int(j) for j in np.flatnonzero(row)) for row in A]
                clique, nodes = maximum_clique(adj)
                graphs += 1
            s = len(clique) + 1
            score = Fraction(s * s * (mu * det - b) ** n, 4**n * det ** (n + 1))
            if score > best:
                best = score
                bestrow = {
                    "K": K.tolist(),
                    "determinant": det,
                    "diameter": str(Fraction(b, det)),
                    "cluster_size": s,
                    "cluster": [[0] * k] + T[ix[clique]].tolist(),
                    "center_density_squared": str(score),
                    "center_density": math.sqrt(float(score)),
                }
    row = {
        "ambient_dimension": N,
        "ambient_minimum": mu,
        "dimension": n,
        "section_rank": k,
        "diameter_limit": str(limit),
        "signed_permutation_classes": all_grams,
        "relaxed_Gram_types": len(types),
        "graphs_solved": graphs,
        "best": bestrow,
        "seconds": time.time() - started,
    }
    out.append(row)
    print(json.dumps(row), flush=True)
    (P / "global_anticodes.json").write_text(json.dumps(out, indent=2) + "\n")
