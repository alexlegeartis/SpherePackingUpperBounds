"""Exact bounded-diameter dual-lattice cluster search in section rank <= 4."""

from pathlib import Path
import itertools, json, math
import numpy as np
from flint import fmpz_mat

P = Path(__file__).resolve().parent


def maximum_clique(adj, max_nodes=None):
    best = []
    nodes = 0

    # Greedy coloring upper bounds, followed by exact branch and bound.
    def colors(Pset):
        ordered = []
        bounds = []
        color = 0
        while Pset:
            color += 1
            available = Pset
            while available:
                bit = available & -available
                v = bit.bit_length() - 1
                ordered.append(v)
                bounds.append(color)
                Pset ^= bit
                available ^= bit
                available &= ~adj[v]
        return ordered, bounds

    def expand(chosen, remaining):
        nonlocal best, nodes
        nodes += 1
        order, bounds = colors(remaining)
        if max_nodes is not None and nodes > max_nodes:
            raise RuntimeError("clique node cap")
        for j in range(len(order) - 1, -1, -1):
            if len(chosen) + bounds[j] <= len(best):
                return
            v = order[j]
            nextset = remaining & adj[v]
            if nextset:
                expand(chosen + [v], nextset)
            elif len(chosen) + 1 > len(best):
                best = chosen + [v]
            remaining &= ~(1 << v)

    expand([], (1 << len(adj)) - 1)
    return best, nodes


if __name__ == "__main__":
    cases = [
        ("44_C5", np.load(P / "antipode44_0.npz")["K"], 44, 509.61898264, 45),
        (
            "45_A3",
            np.full((3, 3), -2, dtype=np.int64) + 8 * np.eye(3, dtype=np.int64),
            45,
            1243.45967883,
            40,
        ),
        ("46_A2", np.array([[6, -3], [-3, 6]]), 46, 2719.93645952, 32),
    ]
    out = []
    for name, K, n, record, bmax in cases:
        inv = fmpz_mat(K.tolist()).inv()
        den = math.lcm(*(int(x.denominator) for row in inv.tolist() for x in row))
        M = np.array(
            [[int(x * den) for x in row] for row in inv.tolist()], dtype=np.int64
        )
        assert np.array_equal(K @ M, den * np.eye(len(K), dtype=np.int64))
        det = int(fmpz_mat(K.tolist()).det())
        maxbeta = bmax / 55 if n == 44 else 0.8
        bounds = [int(math.sqrt(maxbeta * K[i, i])) for i in range(len(K))]
        V = np.array(
            list(itertools.product(*(range(-b, b + 1) for b in bounds))), dtype=np.int64
        )
        norm = np.einsum("vi,ij,vj->v", V, M, V)
        V = V[(norm > 0) & (norm <= maxbeta * den + 1e-7)]
        norm = np.einsum("vi,ij,vj->v", V, M, V)
        D = norm[:, None] + norm[None, :] - 2 * V @ M @ V.T
        for b in sorted(set(map(int, norm))):
            ix = np.flatnonzero(norm <= b)
            A = D[np.ix_(ix, ix)] <= b
            np.fill_diagonal(A, False)
            adj = [sum(1 << int(j) for j in np.flatnonzero(row)) for row in A]
            clique, nodes = maximum_clique(adj)
            S = np.r_[np.zeros((1, len(K)), dtype=np.int64), V[ix[clique]]]
            s = len(S)
            beta = b / den
            density = s * ((6 - beta) / 4) ** (n / 2) / math.sqrt(det)
            row = {
                "name": name,
                "dimension": n,
                "denominator": den,
                "diameter_numerator": b,
                "vertices": len(ix),
                "cluster_size": s,
                "clique_search_nodes": nodes,
                "center_density": density,
                "record_ratio": density / record,
                "cluster": S.tolist(),
            }
            out.append(row)
            print(json.dumps(row), flush=True)
        (P / "dual_anticodes.json").write_text(json.dumps(out, indent=2) + "\n")
