#!/usr/bin/env python3
"""Randomized relaxed Gram screen for non-minimum-basis rank-5 sections.

This is a discovery experiment only. A Gram form is not an embedded section
of P48p until an explicit primitive embedding is supplied and checked.
"""
from itertools import product
import time
import numpy as np
import sympy as sp

SEED = 20261116
SAMPLES = 50_000
rng = np.random.default_rng(SEED)
k = 5
v1 = np.array([v for v in product(range(-1, 2), repeat=k) if any(v)], dtype=np.int64)
v2 = np.array([v for v in product(range(-2, 3), repeat=k) if any(v)], dtype=np.int64)

def exact_minimum(K):
    M = sp.Matrix(K.tolist())
    det = int(M.det())
    inv = M.inv()
    # If v^T K v < 6, then |v_i|^2 < 6 (K^-1)_ii.
    bounds = [int(sp.ceiling(sp.sqrt(6 * inv[i, i]))) + 1 for i in range(k)]
    best, witness, checked = 10**9, None, 0
    for v in product(*[range(-b, b + 1) for b in bounds]):
        if not any(v):
            continue
        q = sum(int(v[i]) * int(K[i, j]) * int(v[j])
                for i in range(k) for j in range(k))
        checked += 1
        if q < best:
            best, witness = q, v
    return det, best, witness, bounds, checked

def main():
    start = time.time()
    n_pd = n1 = n2 = 0
    candidates = []
    for offset in range(0, SAMPLES, 250):
        m = min(250, SAMPLES - offset)
        Ks = np.zeros((m, k, k), dtype=np.int64)
        diag = rng.choice([6, 8, 10, 12, 14, 16, 18], size=(m, k),
                          p=[.08, .25, .25, .20, .12, .07, .03])
        for i in range(k):
            Ks[:, i, i] = diag[:, i]
        ii, jj = np.tril_indices(k, -1)
        vals = rng.integers(-7, 8, size=(m, len(ii)))
        Ks[:, ii, jj] = vals
        Ks[:, jj, ii] = vals

        positive = np.linalg.eigvalsh(Ks.astype(float))[:, 0] > 1e-6
        n_pd += int(positive.sum())
        Ks = Ks[positive]
        if not len(Ks):
            continue
        pass1 = np.einsum("vi,bij,vj->bv", v1, Ks, v1,
                          optimize=True).min(axis=1) >= 6
        n1 += int(pass1.sum())
        Ks = Ks[pass1]
        if not len(Ks):
            continue
        pass2 = np.einsum("vi,bij,vj->bv", v2, Ks, v2,
                          optimize=True).min(axis=1) >= 6
        n2 += int(pass2.sum())
        for K in Ks[pass2]:
            d = int(sp.Matrix(K.tolist()).det())
            if 0 < d < 1800:
                candidates.append((d, K.copy()))

    candidates.sort(key=lambda x: x[0])
    print(f"seed={SEED} samples={SAMPLES}")
    print(f"positive_definite={n_pd} pass_box_1={n1} pass_box_2={n2}")
    print(f"relaxed_candidates_det_lt_1800={len(candidates)}")
    valid = 0
    for d0, K in candidates:
        d, mu, witness, bounds, checked = exact_minimum(K)
        print(f"candidate det={d} exact_min={mu} witness={witness} "
              f"box={bounds} checked={checked}")
        if mu >= 6:
            valid += 1
            print("  PASS exact minimum; embedding into P48p remains unproved")
        else:
            print("  REJECT exact short vector")
        print(K.tolist())
    print(f"exact_minimum_passes={valid}")
    print(f"elapsed_seconds={time.time() - start:.3f}")
    print("SCOPE: abstract rank-5 Gram forms only; no ambient embedding or "
          "packing certificate was found.")

if __name__ == "__main__":
    main()
