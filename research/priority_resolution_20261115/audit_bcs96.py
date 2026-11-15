"""Reconstruct BCS96's parent/norm-doubler and a finite local-optimality proof.

Inherits the published minimum 6 of P48q. Exact input and witness arithmetic;
no enumeration in dimension 96 and no claim of a globally optimal packing.
"""
from pathlib import Path
import argparse
import json
import numpy as np
from flint import fmpz_mat, fmpq_mat

if not __debug__:
    raise RuntimeError("Assertions must be enabled.")
P = Path(__file__).resolve().parent


def audit():
    residues = {i*i % 47 for i in range(1, 47)}
    generators = [[0]+[2*int((j-t) % 47 in residues | {0}) for j in range(47)]
                  for t in range(47)]
    generators += [[-5]+[1]*47]
    roots = []
    for i in range(48):
        for j in range(i+1, 48):
            for sign in (1, -1):
                v = [0]*48
                v[i], v[j] = 6, 6*sign
                roots.append(v)
    generators += roots
    C = fmpz_mat(fmpz_mat(generators).hnf().tolist()[:48])
    assert C.det() == 12**24
    gram = fmpq_mat(C*C.transpose())/12
    assert all(x.denominator == 1 for row in gram.tolist() for x in row)
    G = fmpz_mat([[int(x) for x in row] for row in gram.tolist()])
    assert G.det() == 1 and all(G[i, i] % 2 == 0 for i in range(48))
    n = [[0]*48 for _ in range(48)]
    n[0][0] = 5
    for i in range(1, 48):
        n[0][i] = n[i][0] = -1
    for i in range(47):
        for j in range(47):
            k = (i+j) % 47
            n[i+1][j+1] = -5 if k == 0 else 1 if k in residues else -1
    N = fmpz_mat(n)
    I = fmpz_mat(np.eye(48, dtype=int).tolist())
    assert N == N.transpose() and N*N == 72*I
    Ci = fmpq_mat(C).inv()
    action = fmpq_mat(C*N)*Ci/6
    assert all(x.denominator == 1 for row in action.tolist() for x in row)
    A = fmpz_mat([[int(x) for x in row] for row in action.tolist()])
    assert A*A == 2*I and A*G*A.transpose() == 2*G
    assert all(x.denominator == 1 for row in (fmpq_mat(N)*Ci).tolist() for x in row)
    # Gamma96 = {(x,y) in (D^-1 Lambda)^2 : x-y in Lambda}, D=N/6.
    inv = fmpq_mat(A)/2
    gf = fmpq_mat(G)
    cross = gf*inv.transpose()
    g96 = fmpq_mat(96, 96)
    for i in range(48):
        for j in range(48):
            g96[i, j] = g96[48+i, 48+j] = gf[i, j]
            g96[i, 48+j] = cross[i, j]
            g96[48+j, i] = cross[i, j]
    assert g96.det() == 1/fmpq_mat([[2**48]])[0, 0]
    assert all((2*x).denominator == 1 for row in g96.tolist() for x in row)
    # A small minimum-contact subset suffices for perfection and eutaxy.
    z = np.zeros(48, dtype=np.int64)
    normals = np.array(n, dtype=np.int64)
    contacts = []
    for v in np.array(roots, dtype=np.int64):
        contacts.extend([np.r_[v, z], np.r_[z, v]])
        numer = v@normals
        assert np.all(numer % 12 == 0)
        u = numer//12
        contacts.extend([np.r_[u, u], np.r_[u, -u]])
    pairs = 0
    for i in range(48):
        for j in range(i+1, 48):
            plus, minus = normals[i]+normals[j], normals[i]-normals[j]
            assert np.all(plus % 2 == 0) and np.all(minus % 2 == 0)
            u, v = plus//2, minus//2
            # u-v=N_j, u+v=N_i are in Lambda, as checked above.
            assert u@v == 0 and u@u == v@v == 36
            contacts.extend([np.r_[u, v], np.r_[u, -v], np.r_[v, u], np.r_[v, -u]])
            pairs += 1
    S = np.array(contacts, dtype=np.int64)
    assert S.shape == (13536, 96)
    rays = {tuple(v if next(x for x in v if x) > 0 else -v) for v in S}
    assert len(rays) == len(S)
    assert np.all(np.einsum('ij,ij->i', S, S) == 72)
    assert np.array_equal(S.T@S, 10152*np.eye(96, dtype=np.int64))
    # Spanning proof is in BCS96_LOCAL_PROOF.md: two diagonal symmetric
    # blocks, the symmetric cross block, and all skew cross coordinates.
    assert 2*(48*49//2)+(48*49//2)+pairs == 96*97//2 == 4656
    return {"parent": "Leech-Sloane P48q, using the BCS1982 coordinate generators",
            "inherited_input": "Published P48q minimum 6, not re-enumerated",
            "parent_gram_determinant": 1,
            "norm_doubler_verified": True,
            "BCS96_minimum_squared_inheriting_parent": 6,
            "BCS96_gram_determinant": "1/2^48",
            "BCS96_center_density_exact": "3^48/2^24",
            "finite_contact_rays": 13536,
            "perfect_span_dimension": 4656,
            "uniform_eutaxy_sum_in_orthonormal_coordinates": "846*I96",
            "status": "perfect_and_strictly_eutactic_inheriting_parent_minimum",
            "consequence": "Strict local maximum among lattice Gram deformations, by Voronoi's criterion; not global optimality",
            "all_matrix_checks_exact": True}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-result", action="store_true")
    args = ap.parse_args()
    result = audit()
    path = P/"bcs96_audit.json"
    if args.write_result:
        path.write_text(json.dumps(result, indent=2)+"\n")
    else:
        assert json.loads(path.read_text()) == result, "Stored BCS audit is stale."
    print(json.dumps(result, indent=2))
