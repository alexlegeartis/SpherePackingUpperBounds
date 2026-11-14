"""Floating local metric screens on primitive sections of the 36D record.

Every section basis and inherited minimum vector is checked over the integers.
The metric LP and completeness of the inherited contact shell are numerical.
"""
from pathlib import Path
import itertools, json, math, subprocess, sys
import numpy as np
from fpylll import IntegerMatrix, GSO, LLL

P = Path(__file__).resolve().parent


def bezout(a, b):
    old_r, r, old_s, s, old_t, t = a, b, 1, 0, 0, 1
    while r:
        q = old_r // r
        old_r, r, old_s, s, old_t, t = r, old_r-q*r, s, old_s-q*s, t, old_t-q*t
    if old_r < 0:
        return -old_r, -old_s, -old_t
    return old_r, old_s, old_t


def kernel(h):
    h = h.copy()
    T = np.eye(len(h), dtype=np.int64)
    for i in range(1, len(h)):
        a, b = int(h[0]), int(h[i])
        if not b:
            continue
        g, s, t = bezout(a, b)
        T[:, [0, i]] = T[:, [0, i]] @ np.array([[s, -b//g], [t, a//g]], dtype=np.int64)
        h[0], h[i] = g, 0
    assert h[0] != 0 and not h[1:].any()
    return T[:, 1:].T


def main():
    x = np.load(P/'new36.npz')
    Q = x['Q']
    original_change = np.load(P/'new36_basis_from_original.npz')['U']
    dual = np.load(P/'new36_dual_bkz.npz')
    H = dual['U'] @ original_change.T
    Ad = dual['A'] / float(dual['denom'])
    assert np.max(abs(H @ np.linalg.inv(Q/20) @ H.T-Ad)) < 1e-7
    V = np.loadtxt(P/'new36_vectors.txt', dtype=np.int64)
    assert np.all(np.einsum('vi,ij,vj->v', V, Q, V) == 160)
    choices = [(i,) for i in (0, 6, 10)]
    pairs = sorted(itertools.combinations(range(36), 2), key=lambda ij: np.linalg.det(Ad[np.ix_(ij,ij)]))
    seen = set()
    for ij in pairs:
        d = round(float(np.linalg.det(Ad[np.ix_(ij,ij)])), 7)
        if d in seen:
            continue
        seen.add(d)
        choices.append(ij)
        if len(seen) == 4:
            break
    rows = []
    for selected in choices:
        B = np.eye(36, dtype=np.int64)
        for i in selected:
            B = kernel(B @ H[i]) @ B
        assert not (B @ H[list(selected)].T).any()
        S = B @ Q @ B.T
        b = IntegerMatrix.from_matrix(S.tolist()); u = IntegerMatrix.identity(len(S))
        gso = GSO.Mat(b, U=u, gram=True, flags=GSO.INT_GRAM)
        LLL.Reduction(gso)()
        U = np.array([[int(u[i,j]) for j in range(len(S))] for i in range(len(S))], dtype=np.int64)
        B = U @ B
        S = B @ Q @ B.T
        assert np.array_equal(S, np.array([[int(b[i,j]) for j in range(len(S))] for i in range(len(S))]))
        selected_V = V[np.all(V @ H[list(selected)].T == 0, axis=1)]
        W = np.rint(selected_V @ B.T @ np.linalg.inv(B @ B.T)).astype(np.int64)
        assert np.array_equal(W @ B, selected_V)
        assert np.all(np.einsum('vi,ij,vj->v', W, S, W) == 160)
        prefix = P / ('new36_section' + str(len(S)) + '_' + '_'.join(map(str, selected)))
        np.savez(str(prefix)+'.npz', Q=S, B=B, dual_indices=np.array(selected))
        np.savetxt(str(prefix)+'_vectors.txt', W, fmt='%d')
        delta = 2.**(len(S)/2) / math.exp(np.linalg.slogdet(S/20)[1]/2)
        row = {'dimension':len(S), 'dual_indices':selected, 'center_density':delta,
               'inherited_contacts':len(W)}
        result = subprocess.run([sys.executable, str(P/'eutaxy.py'), str(prefix)], capture_output=True, text=True)
        row['metric_lp_log'] = result.stdout
        row['metric_lp_error'] = result.stderr
        print(json.dumps(row), flush=True)
        rows.append(row)
        (P/'section_metric_screen.json').write_text(json.dumps(rows,indent=2)+'\n')


if __name__ == '__main__':
    main()
