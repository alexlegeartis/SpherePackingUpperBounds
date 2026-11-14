"""Target maximal self-orthogonal neighbors through low-cost hyperplanes."""

from pathlib import Path
import argparse, json, time
from collections import deque
import numpy as np

P = Path(__file__).resolve().parent
from finite import geometry_costs


def rref(M):
    A = M.copy() % 3
    rank = 0
    piv = []
    for j in range(A.shape[1]):
        z = np.flatnonzero(A[rank:, j])
        if not len(z):
            continue
        i = rank + int(z[0])
        A[[rank, i]] = A[[i, rank]]
        A[rank] = A[rank] * int(A[rank, j]) % 3
        for i in range(len(A)):
            if i != rank:
                A[i] = (A[i] - A[i, j] * A[rank]) % 3
        piv.append(j)
        rank += 1
        if rank == len(A):
            break
    return A[:rank], piv


def digits(k):
    return np.array(
        [[q // 3**i % 3 for i in range(k)] for q in range(3**k)], dtype=np.int16
    )


def ft3(C):
    a = C.astype(np.int64).copy()
    b = np.zeros_like(a)
    k = 0
    while 3**k < len(C):
        k += 1
    assert 3**k == len(C)
    for i in range(k):
        A = a.reshape(-1, 3, 3**i)
        B = b.reshape(-1, 3, 3**i)
        a0, a1, a2 = (A[:, j].copy() for j in range(3))
        b0, b1, b2 = (B[:, j].copy() for j in range(3))
        A[:, 0] = a0 + a1 + a2
        B[:, 0] = b0 + b1 + b2
        A[:, 1] = a0 - b1 + b2 - a2
        B[:, 1] = b0 + a1 - b1 - a2
        A[:, 2] = a0 + b1 - a1 - b2
        B[:, 2] = b0 - a1 + a2 - b2
    assert not b.any()
    return a


def extensions(H):
    A, piv = rref(H)
    n = A.shape[1]
    free = [i for i in range(n) if i not in piv]
    N = np.zeros((len(free), n), dtype=np.int16)
    for index, j in enumerate(free):
        N[index, j] = 1
        N[index, piv] = -A[:, j] % 3
    assert not (H @ N.T % 3).any()
    base = A.copy()
    R = []
    for v in N:
        candidate, _ = rref(np.vstack([base, v]))
        if len(candidate) > len(base):
            R.append(v)
            base = candidate
    R = np.array(R)
    D = digits(len(R))
    highest = np.array([next((int(v) for v in reversed(row) if v), 0) for row in D])
    V = D[highest == 1] @ R % 3
    V = V[np.sum(V * V, axis=1) % 3 == 0]
    assert len(V) == (4 if n == 17 else 10)
    return V


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("n", type=int, choices=(17, 18))
    parser.add_argument("--seconds", type=float, default=120)
    parser.add_argument("--hyperplanes", type=int, default=70)
    parser.add_argument("--output-prefix", type=Path)
    parser.add_argument("--seed-file", type=Path)
    args = parser.parse_args()
    n = args.n
    prefix = (args.output_prefix or P / f"work/neighbor{n}").resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    D8 = digits(8)
    D7 = digits(7)
    highest = np.array([next((int(v) for v in reversed(row) if v), 0) for row in D8])
    projective = np.flatnonzero(highest == 1)
    geometric = geometry_costs(n)
    cost = np.where(geometric < 9, 1000, np.where(geometric < 12, 1, 0))
    cost[0] = 0
    encode = 1 << np.arange(n, dtype=np.uint32)
    rng = np.random.default_rng(20261114 + n)
    root = P.parent
    G = np.loadtxt(
        args.seed_file or root / f"results/selforth12_{n}_best.txt",
        skiprows=1,
        dtype=np.int16,
    )
    G, _ = rref(G)
    start = time.time()
    tabu = deque(maxlen=300)
    tabu.append(G.tobytes())
    best = 10**9
    rows = []
    evaluated = 0

    def masks(T):
        return (T != 0).astype(np.uint32) @ encode

    assert G.shape == (8, n) and not (G @ G.T % 3).any()

    def save(H, score, step):
        nonlocal best
        if score >= best:
            return
        assert int(cost[masks(D8 @ H % 3)].sum()) == 2 * score
        best = score
        Path(str(prefix) + "_best.txt").write_text(
            f"{n} 8 {best} {step} 1000\n"
            + "\n".join(" ".join(map(str, row)) for row in H)
            + "\n"
        )
        print("BEST", best, "step", step, "seconds", time.time() - start, flush=True)

    for step in range(10000):
        T = D8 @ G % 3
        C = cost[masks(T)]
        energy = int(C.sum() // 2)
        save(G, energy, step)
        if not best:
            break
        F = ft3(C)
        hyper = (int(C.sum()) + 2 * F) // 3
        assert not np.any((int(C.sum()) + 2 * F) % 3)
        if step == 0:
            for q in rng.choice(projective, 10, replace=False):
                assert int(C[(D8 @ D8[q]) % 3 == 0].sum()) == int(hyper[q])
        order = projective[np.argsort(hyper[projective], kind="stable")]
        selected = list(order[: args.hyperplanes])
        selected += list(rng.choice(order[args.hyperplanes :], 10, replace=False))
        choice = None
        ties = 0
        seen = set()
        for index in selected:
            if not best or time.time() - start > args.seconds:
                break
            h = D8[index]
            pivot = int(np.flatnonzero(h)[0])
            others = [j for j in range(8) if j != pivot]
            K = np.eye(8, dtype=np.int16)[others]
            K[:, pivot] = -h[others] * h[pivot] % 3
            H = K @ G % 3
            T7 = D7 @ H % 3
            inside = int(cost[masks(T7)].sum() // 2)
            for v in extensions(H):
                candidate, _ = rref(np.vstack([H, v]))
                key = candidate.tobytes()
                if key in tabu or key in seen:
                    continue
                seen.add(key)
                score = inside + int(cost[masks((T7 + v) % 3)].sum())
                evaluated += 1
                save(candidate, score, step)
                if choice is None or score < choice[0]:
                    choice = (score, candidate)
                    ties = 1
                elif score == choice[0]:
                    ties += 1
                    if rng.integers(ties) == 0:
                        choice = (score, candidate)
        row = {
            "step": step,
            "energy": energy,
            "minimum_hyperplane_energy": int(hyper[projective].min() // 2),
            "next_energy": None if choice is None else choice[0],
            "evaluated_neighbors": evaluated,
            "seconds": time.time() - start,
        }
        rows.append(row)
        if step % 5 == 0:
            print(row, flush=True)
        if not best or time.time() - start > args.seconds or choice is None:
            break
        G = choice[1]
        tabu.append(G.tobytes())
    Path(str(prefix) + ".json").write_text(
        json.dumps(
            {"best_energy": best, "evaluated_neighbors": evaluated, "iterations": rows},
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
