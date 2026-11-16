"""Discovery of strong 55D projections of MW56; exact witnesses bound the search."""

from pathlib import Path
from fractions import Fraction as F
import numpy as np, json, time, math
from flint import fmpz_mat
from fpylll import IntegerMatrix, GSO, LLL, BKZ, FPLLL
import argparse

P = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser(
    description="Finite projected-coset discovery; output is not a packing certificate."
)
ap.add_argument("--output", type=Path, required=True)
ap.add_argument("--limit", type=int, default=24)
args = ap.parse_args()
if args.output.exists():
    raise SystemExit("Output exists; choose a fresh path.")
if not 1 <= args.limit <= 24:
    raise SystemExit("--limit must be between 1 and 24.")
FPLLL.set_precision(160)
d = np.load(P / "cm58_20261116/certificate.npz", allow_pickle=False)
G = fmpz_mat(d["G"].tolist())
H = 29 * G.inv()
assert all(x.denominator == 1 for row in H.tolist() for x in row)
H = fmpz_mat([[int(x) for x in row] for row in H.tolist()])
roots = d["point_coordinates"]
target = F(5, 2) ** 55 / (2 * 3**24)


def color_upper(adj, neighbors):
    # Greedy proper coloring of the possible-neighbor graph. Every clique
    # can be translated to contain 0, hence this plus1 bounds its size.
    colors = []
    for v in sorted(
        neighbors, key=lambda v: -(adj[v] & sum(1 << j for j in neighbors)).bit_count()
    ):
        for c in range(len(colors)):
            if adj[v] & colors[c] == 0:
                colors[c] |= 1 << v
                break
        else:
            colors.append(1 << v)
    return 1 + len(colors), colors


saved = np.load(P / "mw56_projection55_20261116/certificate.npz", allow_pickle=False)
normals = [(int(m), v) for m, v in zip(saved["normal_numerators"], saved["normals"])][
    : args.limit
]
reports = []
start = time.time()
for case, (m, b) in enumerate(normals):
    bf = fmpz_mat([b.tolist()])
    v0 = bf * H
    assert (bf * H * bf.transpose())[0, 0] == m
    _, S = v0.transpose().hnf(transform=True)
    assert (S * v0.transpose()).tolist() == [[1]] + [[0]] * 55
    U = S.inv().transpose()
    assert all(x.denominator == 1 for row in U.tolist() for x in row)
    B = fmpz_mat([[int(x) for x in row] for row in U.tolist()[1:]])
    labels = B * bf.transpose()
    Q = m * (B * G * B.transpose()) - 29 * (labels * labels.transpose())
    assert Q.det() == m**54 * 29**7
    A = IntegerMatrix.from_matrix([[int(x) for x in row] for row in Q.tolist()])
    V = IntegerMatrix.identity(55)
    M = GSO.Mat(A, U=V, gram=True, float_type="mpfr", update=True)
    lll = LLL.Reduction(M)
    lll()
    BKZ.Reduction(M, lll, BKZ.Param(block_size=32, max_loops=2))()
    V = fmpz_mat([[int(V[i, j]) for j in range(55)] for i in range(55)])
    B = V * B
    Q = V * Q * V.transpose()
    qq = np.array(Q.tolist(), dtype=np.int64)
    bb = np.array(B.tolist(), dtype=np.int64)
    lab = (bb @ b) % m
    best = np.full(m, 10**15, dtype=np.int64)
    witness = [None] * m

    def update(X, N):
        J = X @ lab % m
        for i in np.argsort(N):
            j = int(J[i])
            nv = int(N[i])
            if j and nv < best[j]:
                best[j] = best[-j] = nv
                witness[j] = X[i].copy()
                witness[-j] = -X[i].copy()

    eye = np.eye(55, dtype=np.int64)
    X = np.vstack([eye, -eye])
    N = np.diag(qq).repeat(2).reshape(55, 2).T.flatten()
    assert np.array_equal(N, np.einsum("vi,ij,vj->v", X, qq, X))
    update(X, N)
    seen = {tuple(row) for row in X}
    pool = X.copy()
    for turn in range(7):
        allnew = []
        norms = np.einsum("vi,ij,vj->v", pool, qq, pool)
        dots = pool @ qq
        for j in range(55):
            for sg in [-1, 1]:
                nn = norms + qq[j, j] + 2 * sg * dots[:, j]
                ids = np.flatnonzero((nn > 0) & (nn <= 14 * m))
                if len(ids):
                    xx = pool[ids].copy()
                    xx[:, j] += sg
                    update(xx, nn[ids])
                    for row in xx:
                        key = tuple(row)
                        if key not in seen:
                            seen.add(key)
                            allnew.append(row)
                        if len(allnew) >= 50000:
                            break
                if len(allnew) >= 50000:
                    break
            if len(allnew) >= 50000:
                break
        if not allnew:
            break
        pool = np.array(allnew, dtype=np.int64)
        if len(seen) > 200000:
            break
    # The short-vector walk is discovery only. Fill any missing class with
    # an exact representative using a Bezout lift of the primitive label.
    col, Sb = bf.transpose().hnf(transform=True)
    assert int(col[0, 0]) == 1
    lift = np.array(Sb.tolist()[0], dtype=np.int64)
    V0 = np.array(v0.tolist()[0], dtype=np.int64)
    parent_witness = []
    for j in range(m):
        if j == 0:
            candidates = roots[(roots @ b) == 0]
            assert len(candidates)
            w = candidates[0]
            nn = 10 * m
        elif witness[j] is not None:
            w = witness[j] @ bb
            nn = int(best[j])
        else:
            w = j * lift
            nn = int(
                m
                * (fmpz_mat([w.tolist()]) * G * fmpz_mat([w.tolist()]).transpose())[
                    0, 0
                ]
                - 29 * int(w @ b) ** 2
            )
            best[j] = nn
        # Keep coefficients moderate, without changing the projection.
        w = w - int(round(int(w @ b) / m)) * V0
        assert int(w @ b) % m == j
        wf = fmpz_mat([w.tolist()])
        exact = int(m * (wf * G * wf.transpose())[0, 0] - 29 * int(w @ b) ** 2)
        assert exact == nn and nn > 0
        parent_witness.append(w.tolist())
    best[0] = 10 * m
    # A Cayley graph must be undirected, including when fallback lifts
    # were needed. Negating a witness gives the opposite residue exactly.
    for j in range(1, m):
        if best[j] < best[-j]:
            best[-j] = best[j]
            parent_witness[-j] = [-int(x) for x in parent_witness[j]]
    assert all(best[j] == best[-j] for j in range(1, m))
    # Exhaust all critical radius thresholds using coloring upper bounds.
    bounds = []
    upper_best = F(0)
    unresolved = []
    for cutoff in sorted(set(int(v) for v in best if 0 < v <= 10 * m)):
        allowed = [j for j in range(1, m) if int(best[j]) >= cutoff]
        adj = []
        for i in range(m):
            adj.append(sum(1 << ((i + j) % m) for j in allowed))
        size, colors = color_upper(adj, allowed)
        d2 = F(size * size) * (F(cutoff, 4 * m) ** 55) / (29**7 * m)
        upper_best = max(upper_best, d2)
        bounds.append(
            {
                "cutoff_numerator": cutoff,
                "clique_upper": size,
                "color_classes": [str(x) for x in colors],
            }
        )
        if d2 > target:
            unresolved.append((cutoff, size, float((d2 / target) ** F(1, 2))))
    result = {
        "normal": b.tolist(),
        "normal_norm_numerator": m,
        "normal_norm_denominator": 29,
        "projected_distance_upper_numerators": best.tolist(),
        "parent_witnesses": parent_witness,
        "color_bounds": bounds,
        "density_upper_squared": str(upper_best),
        "printed_55D_target_squared": str(target),
        "exceeds_target_upper": upper_best > target,
        "short_vector_pool": len(seen),
        "unresolved_thresholds": unresolved,
    }
    reports.append(result)
    args.output.write_text(json.dumps(reports, indent=2) + "\n")
    print(
        "CASE",
        case,
        "m",
        m,
        "pool",
        len(seen),
        "distance-range",
        best[1:].min() / m,
        best[1:].max() / m,
        "upper/target",
        math.sqrt(float(upper_best / target)),
        "unresolved",
        unresolved[-5:],
        "seconds",
        time.time() - start,
        flush=True,
    )
