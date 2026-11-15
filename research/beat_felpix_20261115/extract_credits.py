"""Reconstruct the exact Schur range solutions without an optimization solver."""
from pathlib import Path
from fractions import Fraction as Q
from flint import fmpq_mat,fmpq
import json,sys

HERE=Path(__file__).resolve().parent
if not __debug__:raise RuntimeError('assertions must be enabled')
PACKAGE=HERE.parents[1]/'verifications/dim4/refinement_energy20261115'
sys.path.insert(0,str(PACKAGE/'vendor/D4/verify/stage5_kernelbox/src'))
import kformula as kf

records=[]
for d in [5,6,7,8,9]:
    g=kf.GRAM[d]
    matrix=fmpq_mat([[fmpq(str(g.get((i,k),Q(0)))) for k in range(49)] for i in range(49)])
    reduced,rank=matrix.rref()
    indices=[next(k for k in range(49) if reduced[i,k]) for i in range(rank)]
    principal=fmpq_mat([[matrix[i,k] for k in indices] for i in indices])
    restricted=fmpq_mat([[int(i<7)] for i in indices])
    solution=principal.inv()*restricted
    w=fmpq_mat([[solution[indices.index(i),0] if i in indices else 0] for i in range(49)])
    v=fmpq_mat([[int(i<7)] for i in range(49)])
    assert matrix*w==v
    norm=(v.transpose()*w)[0,0]
    assert norm>0
    records.append(dict(degree=d,beta=str(1/norm),solution=[str(w[i,0]) for i in range(49)]))
published=json.loads((PACKAGE/'gram_credits.json').read_text())
assert {'1':records}==published
print('PASS: all five exact Gram credits reconstructed from the original rays')
