"""Exact 86D/87D antipode algebra using the claimed Sakamoto rank-88 parent.

The parent minimum-8 claim is NOT independently verified here. The output
is conditional on that preprint theorem, as labelled in the catalogue.
"""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import re
import numpy as np
import sympy as sp
from flint import fmpz_mat
from sympy.matrices.normalforms import smith_normal_form

if not __debug__:
    raise RuntimeError("Assertions must be enabled.")

P = Path(__file__).resolve().parent


def audit():
    source = P/"sources/dim88a.html"
    body = source.read_text().split('<a NAME="GRAM"', 1)[1]
    body = body.split('</a><br>', 1)[1].split('<a NAME=', 1)[0]
    z = list(map(int, re.findall(r'-?\d+', re.sub('<[^>]*>', ' ', body))))
    assert len(z) == 88**2
    G = sp.Matrix(88, 88, z)
    g = fmpz_mat([[int(x) for x in row] for row in G.tolist()])
    assert G == G.T and g.det() == 1
    assert all(G[i, i] % 2 == 0 for i in range(88))
    assert all(fmpz_mat([row[:j] for row in g.tolist()[:j]]).det() > 0
               for j in range(1, 89))
    rows = []
    for n in (86, 87):
        k = 88-n
        path = P/f"certificates/claimed88_n{n}.npz"
        d = np.load(path, allow_pickle=False)
        assert all(np.issubdtype(d[x].dtype, np.integer) for x in d.files)
        assert sp.Matrix(d["G"].tolist()) == G
        B, K, S = (sp.Matrix(d[x].tolist()) for x in ("B", "K", "S"))
        assert B.shape == (k, 88) and B*G*B.T == K
        assert S == sp.zeros(1, k).col_join(sp.eye(k))
        sf = smith_normal_form(B, domain=sp.ZZ)
        assert all(abs(sf[i, i]) == 1 for i in range(k))
        # Verify surjectivity of integer label lifts and the full period.
        A = fmpz_mat([[int(x) for x in row] for row in (B*G).T.tolist()])
        H, T = A.hnf(transform=True)
        assert T*A == H and abs(T.det()) == 1
        assert H.tolist() == [[int(i == j) for j in range(k)] for i in range(88)]
        N = fmpz_mat(T.tolist()[k:])
        assert (N*g*N.transpose()).det() == int(K.det())
        distances = []
        for i, j in itertools.combinations(range(S.rows), 2):
            v = S.row(i)-S.row(j)
            distances.append((v*K.inv()*v.T)[0])
        beta = max(distances)
        assert beta == {86: sp.Rational(1, 6), 87: sp.Rational(1, 10)}[n]
        assert K.det() == {86: 48, 87: 10}[n]
        d2 = sp.Rational(S.rows**2, K.det())*((8-beta)/4)**n
        rows.append({"dimension": n, "certificate": str(path.relative_to(P)),
                     "normal_gram": K.tolist(), "normal_determinant": str(K.det()),
                     "cosets": S.rows, "squared_diameter": str(beta),
                     "separation_squared_conditional_lower_bound": str(8-beta),
                     "center_density_squared_conditional": str(d2),
                     "center_density_conditional": str(sp.N(sp.sqrt(d2), 40)),
                     "all_matrix_and_cluster_checks_exact": True,
                     "parent_minimum_8_independently_verified": False,
                     "status": "conditional_on_preprint_parent_minimum_theorem",
                     "world_record_confirmed": False,
                     "certificate_sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    return {"source": "Ryotaro Sakamoto, arXiv:2609.15633v2, rank-88 ideal115L0",
            "catalogue_url": "https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/dim88a.html",
            "catalogue_qualification": "The unimodular catalogue labels minimum 8 as claimed by authors",
            "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "parent_minimum_8_independently_verified": False,
            "maximum_dimension": 96, "rows": rows}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-result", action="store_true")
    args = ap.parse_args()
    result = audit()
    serialized = json.loads(json.dumps(result, default=int))
    path = P/"claimed88_comparisons.json"
    if args.write_result:
        path.write_text(json.dumps(serialized, indent=2)+"\n")
    else:
        assert json.loads(path.read_text()) == serialized, "Stored conditional report is stale."
    for row in result["rows"]:
        print(row["dimension"], row["center_density_conditional"], row["status"])
