"""Prepare the five-coset 44D contact search. Discovery, not a certificate."""
import json
import pathlib
import sys
import numpy as np
from fpylll import IntegerMatrix, BKZ

p = pathlib.Path(__file__).resolve().parent / 'dim48mod3'
e = np.load(p / 'antipode44_embedding.npz')
P, w = e['P'], e['w']
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0
rng = np.random.default_rng(seed)
perm = rng.permutation(48)
G = P[np.ix_(perm, perm)]
approx = np.rint(np.linalg.cholesky(G) * 1e8).astype(np.int64)
B = IntegerMatrix.from_matrix(approx.tolist())
U = IntegerMatrix.identity(48)
BKZ.reduction(B, BKZ.Param(block_size=25, max_loops=2), U=U,
              float_type='long double')
U = np.array([[int(U[i,j]) for j in range(48)] for i in range(48)], np.int64)
T = U @ np.eye(48, dtype=np.int64)[perm]
Q = T @ P @ T.T
wf = w @ T.T
np.savez(p / f'periodic44_basis{seed}.npz', Q=Q, T=T, w=wf)
np.savetxt(p / 'periodic44_Q.txt', Q, fmt='%d')
np.savetxt(p / 'periodic44_R.txt', np.linalg.cholesky(Q).T, fmt='%.17g')
np.savetxt(p / 'periodic44_w.txt', wf, fmt='%d')
print(json.dumps({'seed':seed,'diagonal':Q.diagonal().tolist()}),flush=True)
