"""Reproduce a finite minimum-8 discovery pool from the archived L80 Gram.

This is only a seed/closure pool, not an enumeration of the minimum shell.
The saved integer normal certificates do not depend on rerunning this search.
"""
from pathlib import Path
import argparse
import json
import numpy as np
from verify80 import catalogue_gram


def canonical(v):
    v = np.asarray(v, dtype=np.int64)
    v *= np.sign(v[np.arange(len(v)), np.argmax(v != 0, axis=1)])[:, None]
    return np.unique(v, axis=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--rounds", type=int, default=3)
    args = ap.parse_args()
    G = np.array(catalogue_gram(), dtype=np.int64)
    E = np.eye(80, dtype=np.int64)
    seeds = [E[i] for i in range(80) if G[i, i] == 8]
    for i in range(80):
        for j in range(i+1, 80):
            if G[i, j] and G[i, i]+G[j, j]-2*abs(G[i, j]) == 8:
                seeds.append(E[i]-int(np.sign(G[i, j]))*E[j])
    V = canonical(seeds)
    sizes = [len(V)]
    for _ in range(args.rounds):
        ips = V@G@V.T
        i, j = np.where(np.abs(ips) == 4)
        W = V[i]-np.sign(ips[i, j])[:, None]*V[j]
        V = canonical(np.vstack((V, W)))
        sizes.append(len(V))
    assert np.all(np.einsum("vi,vi->v", V@G, V) == 8)
    np.savez_compressed(args.output, G=G, V=V)
    print(json.dumps({"unsigned_pool_sizes": sizes,
                      "all_vector_norms_exact": 8,
                      "full_minimum_shell_enumerated": False}))


if __name__ == "__main__":
    main()
