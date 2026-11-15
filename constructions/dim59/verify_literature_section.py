"""Exact classical kappa13-in-Leech check and derived 59D comparison.

Inherited source inputs: Nebe's Gamma72 is unimodular of minimum 8 and has a
primitive sqrt(2)*Leech section whose orthogonal complement has minimum 8
(arXiv:1008.2862, Section 4.2). The numerical 59D consequence is derived here;
it is not claimed to be a 59D density reported by Nebe.
"""
from pathlib import Path
import argparse
import hashlib
import html
import json
import re
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

if not __debug__:
    raise RuntimeError("Run without -O; exact verification requires assertions.")

P = Path(__file__).resolve().parent


def matrix(name, key):
    text = (P / "comparison_sources" / name).read_text()
    body = text.split('<a NAME="' + key + '"', 1)[1]
    body = body.split('</a><br>', 1)[1].split('<a NAME=', 1)[0]
    values = re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?',
                        html.unescape(re.sub(r'<[^>]*>', ' ', body)))
    m, n = map(int, values[:2])
    assert len(values) == 2 + m*n
    return sp.Matrix(m, n, [sp.Rational(x) for x in values[2:]])


def audit():
    B = matrix("KAPPA13.html", "BASIS")
    H = matrix("KAPPA13.html", "GRAM")
    L = matrix("Leech.html", "BASIS")
    F = matrix("Leech.html", "GRAM")
    assert B.shape == (13, 24) and L.shape == (24, 24)
    # The catalogue coordinate matrices use scale 1/sqrt(8). Their entries
    # are exact integers printed in decimal scientific notation, not rounded
    # irrational triangular bases.
    assert all(x.q == 1 for x in B) and all(x.q == 1 for x in L)
    assert B*B.T/8 == H and L*L.T/8 == F
    assert F == F.T and F.det() == 1 and H.det() == 972
    C = B*L.inv()
    assert all(x.q == 1 for x in C) and C*L == B
    assert C*F*C.T == H
    invariants = smith_normal_form(C, domain=sp.ZZ)
    assert all(abs(invariants[i, i]) == 1 for i in range(13))

    # Primitive embeddings compose. The primitive normal sqrt(2)*kappa13
    # inside Gamma72 has the following determinant. In a unimodular lattice,
    # its primitive perpendicular section has the same Gram determinant.
    normal_det = (2*H).det()
    assert normal_det == 2**13 * 972
    baseline2 = sp.Rational(2**59, normal_det)
    assert baseline2 == sp.Rational(2**44, 3**5)
    ours = json.loads((P / "result.json").read_text())
    assert ours["dimension"] == 59 and ours["all_checks_exact"]
    our2 = sp.Rational(ours["center_density_squared_exact"])
    assert our2 == sp.Rational(3**33*401**59, 2**296*5**2*13**61)
    ratio2 = our2/baseline2
    assert ratio2 > 1
    for lower, upper in [("2.028457", "2.028458")]:
        assert (1+sp.Rational(lower)/100)**2 < ratio2
        assert ratio2 < (1+sp.Rational(upper)/100)**2
    volume = sp.pi**sp.Rational(59, 2)/sp.gamma(sp.Rational(61, 2))
    return {
        "status": "verified_classical_section_derived_from_published_inputs",
        "dimension": 59,
        "name": "Gamma72 perpendicular to sqrt(2)*kappa13",
        "source_ids": ["N12", "K13"],
        "is_reported_59D_result": False,
        "provenance": "The 59D value is derived here from published Gamma72 structure and the classical catalogue kappa13 embedding. It is not quoted as a result of either source.",
        "inherited_theorems": {
            "Gamma72": "Even unimodular, minimum 8; primitive sqrt(2)*Leech section K1 with perpendicular K2 of minimum 8",
            "locator": "Nebe, arXiv:1008.2862, Theorem 2.3 and Section 4.2, pp. 8–9",
            "url": "https://arxiv.org/pdf/1008.2862",
            "why_minimum_is_exact": "The 59D perpendicular contains K2, whose minimum is 8; Gamma72 has no shorter vectors.",
        },
        "kappa13_gram": H.tolist(),
        "kappa13_in_catalogue_Leech_basis": C.tolist(),
        "kappa13_embedding_primitive": True,
        "kappa13_determinant": "972",
        "normal_determinant": str(normal_det),
        "section_determinant": str(normal_det),
        "section_minimum_squared_exact_inheriting_published_theorems": 8,
        "center_density_exact": "2^22/(9*sqrt(3))",
        "center_density_squared_exact": str(baseline2),
        "center_density": str(sp.N(sp.sqrt(baseline2), 40)),
        "packing_fraction": str(sp.N(sp.sqrt(baseline2)*volume, 40)),
        "ratio_squared_exact": str(ratio2),
        "increase_percent": str(sp.N(100*(sp.sqrt(ratio2)-1), 35)),
        "gain_percent_interval_exact_checks": ["2.028457", "2.028458"],
        "all_matrix_and_gain_checks_exact": True,
        "literature_status_report": "../../research/covering_and_deep_cosets_20261115/PRIORITY.md",
        "source_urls": {
            "KAPPA13.html": "https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/KAPPA13.html",
            "Leech.html": "https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/Leech.html",
        },
        "archived_source_normalization": "Trailing whitespace removed and one final newline retained; numerical source data unchanged.",
        "source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in sorted((P / "comparison_sources").glob("*.html"))},
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-result", action="store_true")
    args = parser.parse_args()
    result = audit()
    output = P / "literature_section_comparison.json"
    if args.write_result:
        output.write_text(json.dumps(result, indent=2, default=int) + "\n")
    else:
        assert json.loads(output.read_text()) == result, "Comparison report is stale."
    print(json.dumps({k: result[k] for k in (
        "status", "center_density_exact", "center_density", "packing_fraction",
        "increase_percent", "is_reported_59D_result", "literature_status_report",
    )}, indent=2))
