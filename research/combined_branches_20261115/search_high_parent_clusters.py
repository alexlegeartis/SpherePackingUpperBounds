"""Combine the exact anticode search with the branches' Ne64/Gamma72 parents.

Scores above a reference would be abstract section candidates until an actual
primitive embedding is supplied. No world-record claim follows from this file.
"""

from pathlib import Path
from fractions import Fraction
import itertools, json, math, sys, time
import numpy as np
from flint import fmpz_mat

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(
    0, str(ROOT / "research/lower_bound_search_20261114/broad_dimensions_followup")
)
from dual_anticodes import maximum_clique


def search(N, mu, normal_mu, scale, covol2, k, reference):
    # The physical normal Gram is K/scale. The dual diameter is scale*h*K^-1*h.
    started = time.time()
    ij = np.triu_indices(k, 1)
    seen = set()
    types = []
    signs = [np.array((1,) + s) for s in itertools.product((-1, 1), repeat=k - 1)]
    perms = list(itertools.permutations(range(k)))
    small = np.array(
        [v for v in itertools.product((-1, 0, 1), repeat=k) if any(v)], dtype=np.int64
    )
    for key in itertools.product(
        range(-normal_mu // 2, normal_mu // 2 + 1), repeat=k * (k - 1) // 2
    ):
        if key in seen:
            continue
        K = np.eye(k, dtype=np.int64) * normal_mu
        K[ij] = key
        K[ij[1], ij[0]] = key
        for perm in perms:
            B = K[np.ix_(perm, perm)]
            for s in signs:
                seen.add(tuple(map(int, (s[:, None] * B * s[None, :])[ij])))
        if np.any(np.einsum("vi,ij,vj->v", small, K, small) < normal_mu):
            continue
        if any(fmpz_mat(K[:j, :j].tolist()).det() <= 0 for j in range(1, k + 1)):
            continue
        types.append(K)
    bound = math.isqrt(normal_mu // scale)  # beta <= 1.
    V = np.array(
        [v for v in itertools.product(range(-bound, bound + 1), repeat=k) if any(v)],
        dtype=np.int64,
    )
    best = reference
    improvements = []
    unknown = []
    invalid = []
    graphs = 0
    n = N - k
    for typeid, K in enumerate(types):
        F = fmpz_mat(K.tolist())
        det = int(F.det())
        M = np.array((F.inv() * det).tolist(), dtype=np.int64)
        norm = np.einsum("vi,ij,vj->v", V, M, V)
        keep = scale * norm <= det
        T = V[keep]
        norm = norm[keep]
        order = np.argsort(norm)
        T = T[order]
        norm = norm[order]
        checked = False
        for b in [0] + list(map(int, np.unique(norm))):
            j = int(np.searchsorted(norm, b, side="right"))
            upper = Fraction(
                (j // 2 + 1) ** 2 * scale**k * (mu * det - scale * b) ** n,
                covol2 * 4**n * det ** (n + 1),
            )
            if upper <= best:
                continue
            X = T[:j]
            D = norm[:j, None] + norm[None, :j] - 2 * X @ M @ X.T
            A = D <= b
            np.fill_diagonal(A, False)
            adj = [sum(1 << int(v) for v in np.flatnonzero(row)) for row in A]
            try:
                clique, nodes = maximum_clique(adj, 250000)
            except RuntimeError:
                unknown.append([typeid, str(Fraction(scale * b, det))])
                continue
            graphs += 1
            s = len(clique) + 1
            score = Fraction(
                s * s * scale**k * (mu * det - scale * b) ** n,
                covol2 * 4**n * det ** (n + 1),
            )
            if score > best:
                if not checked:
                    bounds = [
                        math.isqrt((normal_mu - 2) * int(M[i, i]) // det)
                        for i in range(k)
                    ]
                    box = np.array(
                        list(itertools.product(*(range(-v, v + 1) for v in bounds))),
                        dtype=np.int64,
                    )
                    values = np.einsum("vi,ij,vj->v", box, K, box)
                    bad = np.flatnonzero((values > 0) & (values < normal_mu))
                    checked = True
                    if len(bad):
                        invalid.append(
                            {
                                "type": typeid,
                                "vector": box[bad[0]].tolist(),
                                "norm": int(values[bad[0]]),
                            }
                        )
                        break
                best = score
                row = {
                    "normal_Gram_numerator": K.tolist(),
                    "normal_Gram_denominator": scale,
                    "diameter": str(Fraction(scale * b, det)),
                    "cluster": [[0] * k] + X[clique].tolist(),
                    "density_squared": str(score),
                    "ratio_to_reference": math.sqrt(float(score / reference)),
                }
                improvements.append(row)
                print("ABSTRACT CANDIDATE", n, json.dumps(row), flush=True)
    return {
        "dimension": n,
        "parent_dimension": N,
        "parent_minimum": mu,
        "normal_minimum_numerator": normal_mu,
        "normal_scale": scale,
        "covolume_squared": str(covol2),
        "Gram_types": len(types),
        "graphs_solved": graphs,
        "diameter_limit": "1",
        "unknown_graphs": unknown,
        "invalid_Gram_candidates": invalid,
        "best_squared": str(best),
        "reference_squared": str(reference),
        "improvements": improvements,
        "seconds": time.time() - started,
        "world_record_established": False,
    }


if __name__ == "__main__":
    cases = [
        (72, 8, 8, 1, 1, 3, Fraction(4, 75) * Fraction(39, 20) ** 69),
        (72, 8, 8, 1, 1, 4, Fraction(1, 80) * Fraction(39, 20) ** 68),
        (64, 12, 12, 3, 3**32, 2, Fraction(1, 4 * 3**31) * Fraction(35, 12) ** 62),
        (64, 12, 12, 3, 3**32, 3, Fraction(1, 64 * 3**29) * Fraction(93, 32) ** 61),
    ]
    results = []
    for case in cases:
        row = search(*case)
        results.append(row)
        print(json.dumps(row), flush=True)
        (HERE / "high_parent_cluster_results.json").write_text(
            json.dumps(results, indent=2) + "\n"
        )
