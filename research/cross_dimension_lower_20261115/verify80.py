"""Exact audit of four explicit antipode packings from published L80.

The Bachoc--Nebe minimum-8 theorem is inherited, not re-enumerated.
All matrices, saturation, period determinants and diameter inequalities
are checked exactly. FLINT handles large integer matrices; SymPy handles
normal Smith forms and rational calculations. See README.md for proof.
"""
from pathlib import Path
import argparse
import hashlib
import html
import itertools
import json
import re
import numpy as np
import sympy as sp
from flint import fmpz_mat
from sympy.matrices.normalforms import smith_normal_form

if not __debug__:
    raise RuntimeError("Run without -O.")

P = Path(__file__).resolve().parent


def catalogue_gram():
    source = (P/"sources/L80.html").read_text()
    body = source.split('<a NAME="GRAM"', 1)[1].split('</a><br>', 1)[1]
    body = body.split('<a NAME=', 1)[0]
    values = list(map(int, re.findall(r'-?\d+',
                                    html.unescape(re.sub('<[^>]*>', ' ', body)))))
    n, form = values[:2]
    assert n == 80 and form == 0 and len(values) == 2+n*(n+1)//2
    g = [[0]*n for _ in range(n)]
    pos = 2
    for i in range(n):
        for j in range(i+1):
            g[i][j] = g[j][i] = values[pos]
            pos += 1
    return g


def audit():
    g = catalogue_gram()
    fg = fmpz_mat(g)
    assert fg.det() == 1
    assert all(g[i][i] % 2 == 0 for i in range(80))
    # Exact Sylvester criterion. No floating Cholesky tolerance is a proof input.
    assert all(fmpz_mat([row[:j] for row in g[:j]]).det() > 0
               for j in range(1, 81))
    G = sp.Matrix(g)
    ceiling2 = {79: sp.Integer(2)**76, 78: sp.Rational(2**74, 3),
                77: sp.Integer(2)**69, 76: sp.Integer(2)**66}
    root_dets = {76: 4, 77: 4, 78: 3, 79: 2}
    rows = []
    for n in range(79, 75, -1):
        k = 80-n
        path = P/f"certificates/L80_n{n}.npz"
        data = np.load(path, allow_pickle=False)
        assert set(data.files) == {"G", "B", "K", "S"}
        assert all(np.issubdtype(data[key].dtype, np.integer) for key in data.files)
        assert data["G"].tolist() == g
        B, K, S = (sp.Matrix(data[key].tolist()) for key in ("B", "K", "S"))
        assert B.shape == (k, 80) and K.shape == (k, k)
        assert B*G*B.T == K
        normal_det = K.det()
        smith = smith_normal_form(B, domain=sp.ZZ)
        assert all(abs(smith[i, i]) == 1 for i in range(k))

        # Complete kernel and right inverse, computed independently of the
        # Smith saturation test. HNF: T*(B*G)^T = [I_k; 0].
        A = fmpz_mat([[int(x) for x in row] for row in (B*G).T.tolist()])
        H, T = A.hnf(transform=True)
        assert T*A == H and abs(T.det()) == 1
        assert H.tolist() == [[int(i == j) for j in range(k)] for i in range(80)]
        N = fmpz_mat(T.tolist()[k:])
        Q = N*fg*N.transpose()
        assert Q.det() == int(normal_det)
        assert N*A == fmpz_mat(80-k, k)
        R = sp.Matrix(T.tolist()[:k])
        assert R*G*B.T == sp.eye(k)

        inverse = K.inv()
        assert S.shape == (k+1, k)
        assert len(set(tuple(S.row(i)) for i in range(S.rows))) == S.rows
        distances = []
        for i, j in itertools.combinations(range(S.rows), 2):
            v = S.row(i)-S.row(j)
            distances.append((v*inverse*v.T)[0])
        beta = max(distances)
        assert beta < 8
        # For every label, R supplies an integer ambient lift. Subtraction of
        # its normal projection gives the explicit packing shift in the period
        # space. No nearest-point or numerical embedding assumption remains.
        lifts = S*R
        assert lifts*G*B.T == S
        shifts = lifts-S*inverse*B
        assert shifts*G*B.T == sp.zeros(S.rows, k)
        lower_min = 8-beta
        d2 = sp.Rational(S.rows**2, normal_det)*(lower_min/4)**n
        ratio2 = d2/ceiling2[n]
        assert ratio2 > 1
        r = n-72
        product2 = sp.Rational(2**72, 2**r*root_dets[n])
        rows.append({
            "dimension": n, "certificate": str(path.relative_to(P)),
            "normal_gram": K.tolist(), "normal_rank": k,
            "normal_gram_determinant_exact": str(normal_det),
            "normal_primitive": True,
            "period_gram_determinant_exact": str(Q.det()),
            "cosets": S.rows, "dual_squared_diameter_exact": str(beta),
            "packing_minimum_squared_lower_bound_exact": str(lower_min),
            "minimum_attainment_established": False,
            "center_density_squared_exact": str(d2),
            "center_density": str(sp.N(sp.sqrt(d2), 40)),
            "packing_fraction": str(sp.N(sp.sqrt(d2)*sp.pi**sp.Rational(n, 2)/sp.gamma(sp.Rational(n+2, 2)), 40)),
            "bare_section_comparison_scope": "Ceiling only for codimension-k bare-section guarantees that retain the L80 parent minimum 8; not a ceiling for every actual section packing or a published record",
            "bare_section_ceiling_center_density_squared_exact": str(ceiling2[n]),
            "ratio_to_bare_section_ceiling_squared_exact": str(ratio2),
            "increase_over_bare_section_ceiling_percent": str(sp.N(100*(sp.sqrt(ratio2)-1), 30)),
            "product_comparison_method": "Gamma72 orthogonal product with minimum-2 D4/D5/E6/E7, rescaled to minimum 8",
            "product_center_density_squared_exact": str(product2),
            "ratio_to_product_squared_exact": str(d2/product2),
            "published_record_improvement_established": False,
            "all_matrix_and_diameter_checks_exact": True,
            "source_sha256": hashlib.sha256((P/"sources/L80.html").read_bytes()).hexdigest(),
            "certificate_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        })
    return {
        "status": "exact_explicit_antipode_constructions",
        "inherited_input": "Bachoc--Nebe L80: even unimodular, minimum 8; minimum theorem not re-enumerated",
        "source_url": "https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/L80.html",
        "source_paper": "C. Bachoc and G. Nebe, Extremal lattices of minimum 8 related to the Mathieu group M22, J. reine angew. Math. 494 (1998), 155-171",
        "source_archive_normalization": "Trailing whitespace removed and a final newline retained; numerical data unchanged",
        "rows": rows,
        "current_record_identified": False,
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-result", action="store_true")
    args = ap.parse_args()
    result = audit()
    path = P/"l80_results.json"
    if args.write_result:
        path.write_text(json.dumps(result, indent=2, default=int)+"\n")
    else:
        assert json.loads(path.read_text()) == json.loads(json.dumps(result, default=int)), "Stored L80 report is stale."
    for row in result["rows"]:
        print(f'{row["dimension"]}D: center density={row["center_density"]}, '
              f'gain over METHOD ceiling={row["increase_over_bare_section_ceiling_percent"]}%')
