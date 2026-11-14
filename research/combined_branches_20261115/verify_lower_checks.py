"""Exact fixed-cluster maxima and integer upper-distance witness checks."""

from pathlib import Path
import itertools, json, sys
import numpy as np
from flint import fmpz_mat

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(
    0, str(ROOT / "research/lower_bound_search_20261114/broad_dimensions_followup")
)
from dual_anticodes import maximum_clique

parent = ROOT / "research/novel_routes_20261113/higher_parent"
rows = []
for k in (6, 7, 8):
    d = np.load(parent / f"gamma72_root_rank{k}.npz")
    K = d["K"]
    B = d["B"]
    G = d["G"]
    assert np.array_equal(B @ G @ B.T, K)
    F = fmpz_mat(K.tolist())
    det = int(F.det())
    M = np.array((F.inv() * det).tolist(), dtype=np.int64)
    assert np.array_equal(K @ M, det * np.eye(k, dtype=np.int64))
    V = np.array(list(itertools.product(range(-2, 3), repeat=k)), dtype=np.int64)
    norm = np.einsum("vi,ij,vj->v", V, M, V)
    keep = (norm > 0) & (2 * norm <= det)
    V = V[keep]
    norm = norm[keep]
    distances = norm[:, None] + norm[None, :] - 2 * V @ M @ V.T
    A = 2 * distances <= det
    np.fill_diagonal(A, False)
    adj = [sum(1 << int(j) for j in np.flatnonzero(row)) for row in A]
    clique, nodes = maximum_clique(adj, 1000000)
    assert len(clique) + 1 == 9
    rows.append(
        {
            "dimension": 72 - k,
            "normal_rank": k,
            "diameter": "1/2",
            "eligible_nonzero_points": len(V),
            "maximum_cluster_size": 9,
            "clique_search_nodes": nodes,
        }
    )
d = np.load(HERE / "ne64_order2_screen.npz")
source = np.load(parent / "ne64_primal_best.npz")
Q = source["R"]
T = source["T"]
assert np.array_equal(T @ source["G"] @ T.T, Q)
assert abs(fmpz_mat(T.tolist()).det()) == 1
Y = d["witnesses"].astype(np.int64)
assert np.all((Y - d["residues"]) % 2 == 0)
assert 64**2 * int(abs(Q).max()) * int(abs(Y).max()) ** 2 < 2**63
norm = np.einsum("vi,ij,vj->v", Y, Q, Y)
assert np.array_equal(norm, d["norms"])
rejected = int((norm <= 34).sum())
unknown = int((norm > 34).sum())
assert (rejected, unknown) == (12, 88)
report = {
    "all_passed": True,
    "fixed_cluster_checks": rows,
    "Ne64_exact_rejections": rejected,
    "Ne64_capped_inconclusive": unknown,
    "scope": "Integer contact/upper-distance checks and exact finite clique maxima. No new lower-bound record or complete distance certificate.",
}
(HERE / "LOWER_CHECKS.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report), flush=True)
