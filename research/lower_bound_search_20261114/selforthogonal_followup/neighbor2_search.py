"""Enumerate all maximal self-orthogonal completions of targeted dimension-six cuts."""

from pathlib import Path
from collections import deque
import argparse, json, time
import numpy as np
from neighbor_search import rref, digits, ft3

P = Path(__file__).resolve().parent
from finite import geometry_costs


def kernel_cut(G, h):
    k = len(G)
    p = int(np.flatnonzero(h)[0])
    other = [i for i in range(k) if i != p]
    K = np.eye(k, dtype=np.int16)[other]
    K[:, p] = -h[other] * h[p] % 3
    return rref(K @ G % 3)[0]


def projective(D):
    return np.array([next((int(v) for v in reversed(row) if v), 0) == 1 for row in D])


def points_planes(H):
    A, piv = rref(H)
    n = A.shape[1]
    free = [i for i in range(n) if i not in piv]
    N = np.zeros((len(free), n), dtype=np.int16)
    for index, j in enumerate(free):
        N[index, j] = 1
        N[index, piv] = -A[:, j] % 3
    base = A.copy()
    R = []
    for v in N:
        candidate, _ = rref(np.vstack([base, v]))
        if len(candidate) > len(base):
            R.append(v)
            base = candidate
    R = np.array(R)
    d = len(R)
    D = digits(d)
    Q = D[projective(D)]
    V = Q @ R % 3
    keep = np.sum(V * V, axis=1) % 3 == 0
    Q = Q[keep]
    V = V[keep]
    assert len(Q) == (40 if n == 17 else 112)
    power = 3 ** np.arange(d, dtype=np.int64)
    lookup = {int(q @ power): i for i, q in enumerate(Q)}
    orth = V @ V.T % 3
    planes = set()

    def normalize(q):
        return q * next(int(x) for x in reversed(q) if x) % 3

    for i in range(len(Q)):
        for j in range(i + 1, len(Q)):
            if orth[i, j]:
                continue
            a = lookup[int(normalize((Q[i] + Q[j]) % 3) @ power)]
            b = lookup[int(normalize((Q[i] - Q[j]) % 3) @ power)]
            planes.add(tuple(sorted((i, j, a, b))))
    assert len(planes) == (40 if n == 17 else 280)
    return V, np.array(sorted(planes), dtype=np.int16)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("n", type=int, choices=(17, 18))
    parser.add_argument("--seconds", type=float, default=180)
    parser.add_argument("--first", type=int, default=5)
    parser.add_argument("--second", type=int, default=3)
    parser.add_argument("--output-prefix", type=Path)
    parser.add_argument("--seed-file", type=Path)
    args = parser.parse_args()
    n = args.n
    prefix = (args.output_prefix or P / f"work/neighbor2_{n}").resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    Ds = {k: digits(k) for k in (6, 7, 8)}
    reps = {k: np.flatnonzero(projective(Ds[k])) for k in (7, 8)}
    d = geometry_costs(n)
    cost = np.where(d < 9, 1000, np.where(d < 12, 1, 0))
    cost[0] = 0
    encode = 1 << np.arange(n, dtype=np.uint32)
    root = P.parent
    G = rref(
        np.loadtxt(
            args.seed_file or root / f"results/selforth12_{n}_best.txt",
            skiprows=1,
            dtype=np.int16,
        )
    )[0]
    rng = np.random.default_rng(20261130 + n)
    start = time.time()
    tabu = deque([G.tobytes()], maxlen=400)
    rows = []
    evaluated = 0
    best = 10**9

    def masks(T):
        return (T != 0).astype(np.uint32) @ encode

    assert G.shape == (8, n) and not (G @ G.T % 3).any()

    def save(H, score, step):
        nonlocal best
        if score >= best:
            return
        assert int(cost[masks(Ds[8] @ H % 3)].sum()) == 2 * score
        best = score
        Path(str(prefix) + "_best.txt").write_text(
            f"{n} 8 {best} {step} 1000\n"
            + "\n".join(" ".join(map(str, row)) for row in H)
            + "\n"
        )
        print("BEST", best, "step", step, "seconds", time.time() - start, flush=True)

    for step in range(10000):
        T = Ds[8] @ G % 3
        C = cost[masks(T)]
        energy = int(C.sum() // 2)
        save(G, energy, step)
        F = ft3(C)
        scores = (int(C.sum()) + 2 * F) // 3
        order = reps[8][np.argsort(scores[reps[8]], kind="stable")]
        h1s = list(order[: args.first]) + list(
            rng.choice(order[args.first :], 2, replace=False)
        )
        cuts = {}
        for h1 in h1s:
            H7 = kernel_cut(G, Ds[8][h1])
            T7 = Ds[7] @ H7 % 3
            C7 = cost[masks(T7)]
            F7 = ft3(C7)
            s7 = (int(C7.sum()) + 2 * F7) // 3
            o7 = reps[7][np.argsort(s7[reps[7]], kind="stable")]
            h2s = list(o7[: args.second]) + [int(rng.choice(o7[args.second :]))]
            for h2 in h2s:
                H = kernel_cut(H7, Ds[7][h2])
                cuts[H.tobytes()] = H
        choice = None
        ties = 0
        mininside = 10**9
        zeropoints = 0
        for H in cuts.values():
            if time.time() - start > args.seconds:
                break
            T6 = Ds[6] @ H % 3
            inside = int(cost[masks(T6)].sum() // 2)
            mininside = min(mininside, inside)
            V, planes = points_planes(H)
            weights = np.array(
                [int(cost[masks((T6 + v) % 3)].sum()) for v in V], dtype=np.int64
            )
            zeropoints = max(zeropoints, int((weights == 0).sum()))
            energies = inside + weights[planes].sum(axis=1)
            evaluated += len(planes)
            if step == 0:
                for ix in (0, len(planes) - 1):
                    a, b = planes[ix, :2]
                    check = np.vstack([H, V[a], V[b]])
                    assert int(cost[masks(Ds[8] @ check % 3)].sum()) == 2 * int(
                        energies[ix]
                    )
            for ix in np.argsort(energies):
                score = int(energies[ix])
                if choice is not None and score > choice[0]:
                    break
                a, b = planes[ix, :2]
                candidate = rref(np.vstack([H, V[a], V[b]]))[0]
                assert (
                    candidate.shape == (8, n)
                    and not (candidate @ candidate.T % 3).any()
                )
                key = candidate.tobytes()
                if key in tabu:
                    continue
                save(candidate, score, step)
                if choice is None or score < choice[0]:
                    choice = (score, candidate)
                    ties = 1
                else:
                    ties += 1
                    if rng.integers(ties) == 0:
                        choice = (score, candidate)
        row = {
            "step": step,
            "energy": energy,
            "best": best,
            "cuts": len(cuts),
            "minimum_base_energy": mininside,
            "max_zero_cosets": zeropoints,
            "next_energy": None if choice is None else choice[0],
            "completions_evaluated": evaluated,
            "seconds": time.time() - start,
        }
        rows.append(row)
        print(row, flush=True)
        if not best or time.time() - start > args.seconds or choice is None:
            break
        G = choice[1]
        tabu.append(G.tobytes())
    Path(str(prefix) + ".json").write_text(
        json.dumps(
            {
                "best_energy": best,
                "completions_evaluated": evaluated,
                "iterations": rows,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
