"""Construct the classical BCS88 lattice and its primitive block sections.

Exact HNF construction and period checks; minimum 8 inherits BCS1982's
published lattice theorem, whose MDS code inputs are independently checked.
"""
from pathlib import Path
import argparse
import json
import numpy as np
from flint import fmpz_mat, fmpq_mat
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"priority_resolution_20261115"))
from verify_classical_sections import audit as code_audit, mul, power

if not __debug__:
    raise RuntimeError("Assertions must be enabled.")
P = Path(__file__).resolve().parent


def ints(a):
    assert all(x.denominator == 1 for row in a.tolist() for x in row)
    return fmpz_mat([[int(x) for x in row] for row in a.tolist()])


def construct():
    code_audit()
    root = [[int(j == i)-int(j == i+1) for j in range(8)] for i in range(7)]
    root += [[int(j in (6, 7)) for j in range(8)]]
    C8 = fmpz_mat(fmpz_mat([[2*x for x in r] for r in root]+[[1]*8]).hnf().tolist()[:8])
    assert C8.det() == 2**8
    G8 = ints(fmpq_mat(C8*C8.transpose())/4)
    assert G8.det() == 1 and all(G8[i, i] % 2 == 0 for i in range(8))
    D = [[0]*8 for _ in range(8)]
    for i in range(0, 8, 2):
        D[i][i] = D[i][i+1] = D[i+1][i] = 1
        D[i+1][i+1] = -1
    A = ints(fmpq_mat(C8*fmpz_mat(D))*fmpq_mat(C8).inv())
    assert A*A == 2*fmpz_mat(np.eye(8, dtype=int).tolist())
    assert A*G8*A.transpose() == 2*G8
    # Choose a binary basis of Lambda/D Lambda by labels z*A mod 2.
    pivots, basis = [], []
    for i, row in enumerate(A.tolist()):
        v = sum((int(x) % 2) << j for j, x in enumerate(row))
        for p, w in basis:
            if v >> p & 1:
                v ^= w
        if v:
            basis.append((v.bit_length()-1, v))
            pivots.append(i)
    assert len(pivots) == 4
    lifts = [[int(j == i) for j in range(8)] for i in pivots]
    quotient_labels = set()
    for bits in range(16):
        v = fmpz_mat([[sum(((bits >> t) & 1)*lifts[t][i] for t in range(4)) for i in range(8)]])*A
        quotient_labels.add(tuple(int(x) % 2 for x in v.tolist()[0]))
    assert len(quotient_labels) == 16
    raw = (4*fmpz_mat(np.eye(88, dtype=int).tolist())).tolist()
    for level, k in [(1, 10), (2, 8), (3, 4)]:
        transport = 2*A if level == 1 else 2*fmpz_mat(np.eye(8, dtype=int).tolist()) if level == 2 else A
        for j in range(k):
            values = [power(x, j) for x in range(11)]
            for b in (1, 2, 4, 8):
                whole = []
                for value in values:
                    value = mul(b, value)
                    v = [[sum(((value >> t) & 1)*lifts[t][i] for t in range(4)) for i in range(8)]]
                    whole.extend((fmpz_mat(v)*transport).tolist()[0])
                raw.append(whole)
    assert len(raw) == 176
    C = fmpz_mat(fmpz_mat(raw).hnf().tolist()[:88])
    assert C.det() == 2**88
    block = fmpz_mat(88, 88)
    for j in range(11):
        for a in range(8):
            for b in range(8):
                block[8*j+a, 8*j+b] = G8[a, b]
    gram4 = C*block*C.transpose()
    G = fmpq_mat(gram4)/4  # length-scaled by 2 from the min-2 construction
    assert G.det() == 1
    R1 = [[1, -1, 0, 0, 0, 0, 0, 0]]
    R2 = R1+[[0, 1, -1, 0, 0, 0, 0, 0]]
    R3 = R2+[[0, 0, 1, -1, 0, 0, 0, 0]]
    R4 = R3+[[0, 0, 1, 1, 0, 0, 0, 0]]
    R5 = [[int(j == i)-int(j == i+1) for j in range(8)] for i in range(5)]
    D5 = R5[:4]+[[int(j in (3,4)) for j in range(8)]]
    rows = []
    for r, physical in [(5, R5), (5, D5)]:
        roots8 = ints(2*fmpq_mat(physical)*fmpq_mat(C8).inv())
        u = fmpz_mat([row+[0]*80 for row in (roots8*A).tolist()])
        B = ints(fmpq_mat(u*block*C.transpose())/2)
        K = u*block*u.transpose()
        H, T = B.transpose().hnf(transform=True)
        assert H.tolist() == [[int(i == j) for j in range(r)] for i in range(88)]
        assert abs(T.det()) == 1
        N = fmpq_mat(T.tolist()[r:])
        period = N*G*N.transpose()
        assert period.det() == K.det()
        assert all(x == 0 for row in (N*fmpq_mat(B.transpose())).tolist() for x in row)
        w=fmpq_mat([[4*int(i==8) for i in range(88)]])*fmpq_mat(C).inv()
        assert all(x.denominator==1 for row in w.tolist() for x in row)
        assert (w*G*w.transpose())[0,0]==8
        assert all(x==0 for row in (w*fmpq_mat(B.transpose())).tolist() for x in row)
        rows.append({"minimum8_witness_coordinates": [str(x) for x in w.tolist()[0]],"dimension": 88-r, "normal_dual_gram": [[int(x) for x in row] for row in K.tolist()],
                     "normal_primitive": True, "period_gram_determinant": str(period.det())})
    return {"status": "exact_classical_BCS88_construction_and_sections_checked",
            "published_minimum_input": "BCS1982 Theorems 1-2, Eq.(15), E8 example: minimum 2 before length scaling by 2",
            "minimum_squared_after_scaling": 8, "covolume": 1,
            "code_field": "GF(16), x^4+x+1", "code_dimensions": [10, 8, 4],
            "code_evaluation_points": list(range(11)),
            "MDS_maximal_minors_checked": 506,
            "normal_block_dual_minimum": 4, "sections": rows}, gram4, C


if __name__ == '__main__':
    result, gram4, C = construct()
    saved=np.load(Path(__file__).resolve().parents[1]/'priority_resolution_20261115/bcs88_data.npz')
    assert saved['gram_numerator'].tolist()==gram4.tolist()
    assert saved['generator_numerator'].tolist()==C.tolist()
    import mpmath as mp
    mp.mp.dps=80
    for row in result['sections']:
        center=mp.power(2,mp.mpf(83)/2)/mp.sqrt(int(row['period_gram_determinant']))
        density=center*mp.pi**(mp.mpf(83)/2)/mp.gamma(mp.mpf(85)/2)
        row['center_density']=str(center); row['density']=str(density)
        row['gain_percent_current']=str(100*(mp.mpf('2.56962386473426e-18')/density-1))
    Path(__file__).with_name('result.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
