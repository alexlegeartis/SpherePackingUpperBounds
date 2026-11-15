"""Exact independent audit of the Ne64 six-translate 59D construction.

The only inherited minimum theorem is that the archived published Ne64 Gram
has minimum 12. No capped search, floating enumeration or optimizer is trusted.
"""
from pathlib import Path
import argparse
import hashlib
import html
import itertools
import json
import re
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

if not __debug__:
    raise RuntimeError("Exact proof checks require Python assertions; run without -O.")

P = Path(__file__).resolve().parent
ROOT = P.parents[1]
SOURCE = ROOT / "research/novel_routes_20261113/dual64/CQ64b.html"


def read(name):
    return sp.Matrix([[int(x) for x in row.split()]
                      for row in (P / "data" / name).read_text().splitlines()])


def audit():
    text = SOURCE.read_text()
    body = re.search(r'<a NAME="GRAM".*?</a><br>(.*?)<a NAME="DIVISORS"',
                     text, re.S).group(1)
    nums = list(map(int, re.findall(r'-?\d+', html.unescape(re.sub(r'<[^>]*>', ' ', body)))))
    assert nums[:2] == [64, 0] and len(nums) == 2082
    G = sp.zeros(64)
    cursor = 2
    for i in range(64):
        for j in range(i+1):
            G[i, j] = G[j, i] = nums[cursor]
            cursor += 1
    L, D = G.LDLdecomposition(hermitian=False)
    assert L*D*L.T == G and all(D[i,i] > 0 for i in range(64))
    assert sp.prod(D[i,i] for i in range(64)) == 3**32
    assert all(G[i,i] % 2 == 0 for i in range(64))
    Gi = G.inv()

    A, S = read("normal_coefficients.txt"), read("dual_cluster.txt")
    assert A.shape == (5,64) and S == sp.zeros(1,5).col_join(sp.eye(5))
    K = A*Gi*A.T
    expected = 14*sp.eye(5)-2*sp.ones(5)
    for i in range(4): expected[i,i+1] = expected[i+1,i] = -3
    assert 3*K == read("normal_gram_numerator.txt") == expected
    assert K.det() == sp.Rational(33800,81)
    unit = list(map(int, read("unit_minor_columns.txt")))
    assert len(unit) == len(set(unit)) == 5 and all(0 <= c < 64 for c in unit)
    assert abs(A[:,unit].det()) == 1
    snf = smith_normal_form(A, domain=sp.ZZ)
    assert all(abs(snf[i,i]) == 1 for i in range(5))
    free = [c for c in range(64) if c not in unit]
    N, Q = read("section_basis.txt"), read("period_gram.txt")
    assert N.shape == (59,64) and N[:,free] == sp.eye(59)
    assert N*A.T == sp.zeros(59,5)
    # The free identity block and unit minor prove this is the WHOLE integer
    # kernel, so the period contains every ambient point in the section.
    assert Q == N*G*N.T
    period_det = Q.det()
    assert period_det == 3**32*K.det() == 33800*3**28
    lifts = read("ambient_lifts.txt")
    assert lifts.shape == (6,64) and lifts*A.T == S
    denominator = int(read("shift_denominator.txt")[0])
    assert denominator > 0
    shifts = read("shift_numerators.txt") / denominator
    assert shifts.shape == (6,59)
    projected = lifts-S*K.inv()*A*Gi
    assert shifts*N == projected and projected*A.T == sp.zeros(6,5)
    assert shifts[0,:] == sp.zeros(1,59)

    distances = []
    for i,j in itertools.combinations(range(6),2):
        ds = S[i,:]-S[j,:]
        q = (ds*K.inv()*ds.T)[0]
        assert q > 0
        distances.append(q)
        # Equality of two projected cosets would give a nonzero ambient
        # vector wholly in the normal space, of norm q < 12, a contradiction.
        assert q < 12
    beta = max(distances)
    assert beta == sp.Rational(45,104)
    minimum = 12-beta
    assert minimum == sp.Rational(1203,104)

    # A matching ambient norm-12 vector makes the projected minimum exact.
    contact = json.loads((P/"data/contact.json").read_text())
    i,j = contact["i"],contact["j"]
    assert 0 <= i < 6 and 0 <= j < 6 and i != j
    x = sp.Matrix([contact["ambient_coefficients"]])
    ds = S[i,:]-S[j,:]
    assert x.shape == (1,64) and x*A.T == ds and (x*G*x.T)[0] == 12
    assert (ds*K.inv()*ds.T)[0] == beta
    z = x-(lifts[i,:]-lifts[j,:])
    assert z[:,free]*N == z and all(v.q == 1 for v in z)
    c = z[:,free]+shifts[i,:]-shifts[j,:]
    assert (c*Q*c.T)[0] == minimum

    delta2 = 36*(minimum/4)**59/period_det
    assert delta2 == sp.Rational(3**33*401**59,2**296*5**2*13**61)
    baseline2 = sp.Integer(2)**36
    assert delta2 > baseline2
    gain = 100*(sp.sqrt(delta2/baseline2)-1)
    volume = sp.pi**sp.Rational(59,2)/sp.gamma(sp.Rational(61,2))
    result = {
        "status":"verified_periodic_packing", "dimension":59, "cosets":6,
        "normal_rank":5, "normal_gram_numerator":expected.tolist(),
        "normal_gram_denominator":3, "normal_determinant":str(K.det()),
        "cluster_squared_diameter":str(beta), "period_determinant":str(period_det),
        "minimum_squared_exact":str(minimum), "minimum_attained_by_saved_contact":True,
        "center_density_squared_exact":str(delta2),
        "center_density":str(sp.N(sp.sqrt(delta2),35)),
        "packing_fraction":str(sp.N(sp.sqrt(delta2)*volume,35)),
        "comparison_center_density_exact":"262144", "comparison_center_density_squared":"68719476736",
        "comparison_source":"research/novel_routes_20261113/audits/gamma72_leech_section_rank13.npz",
        "comparison_kind":"internal_repository_comparison",
        "comparison_is_published_59D_record":False,
        "comparison_scope":"The Gamma72 rank-13 section was computed in this repository, not identified as a published 59D record. The percentage does not establish an improvement over published SOTA.",
        "ratio_squared_exact":str(delta2/baseline2), "increase_percent":str(sp.N(gain,30)),
        "ambient_minimum_input":12,
        "ambient_minimum_scope":"Published Ne64 theorem is inherited; this audit does not enumerate its 138,458,880 minimum vectors.",
        "ambient_source":"https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/CQ64b.html",
        "ambient_source_sha256":hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "data_sha256":{f.name:hashlib.sha256(f.read_bytes()).hexdigest()for f in sorted((P/"data").iterdir())if f.is_file()},
        "all_checks_exact":True, "world_record_established":False,
    }
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-result",action="store_true")
    args = parser.parse_args()
    result = audit()
    if args.write_result: (P/"result.json").write_text(json.dumps(result,indent=2,default=int)+"\n")
    print(json.dumps({k:result[k]for k in ("status","dimension","cosets","minimum_squared_exact","center_density","packing_fraction","increase_percent","all_checks_exact")},indent=2))
