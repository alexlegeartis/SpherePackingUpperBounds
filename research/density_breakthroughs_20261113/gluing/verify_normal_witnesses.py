"""Exact verification of saved 38D/39D reference-section normal witnesses.

The mother lattice's minimum 6 is a published input, not re-enumerated here.
"""
from pathlib import Path
from fractions import Fraction
import json, math
import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

P = Path(__file__).resolve().parent
rows = []
for rank, expected in [(9, 21870), (10, 54675)]:
    j = np.load(P/f'P48n_best_rank{rank}.npz')
    G, B, K = (sp.Matrix(j[key]) for key in ['G', 'B', 'K'])
    assert G.det() == 1
    assert B*G*B.T == K
    assert K.det() == expected
    S = smith_normal_form(B, domain=sp.ZZ)
    index = math.prod(abs(int(S[i, i])) for i in range(rank))
    assert index == 1
    dimension = 48-rank
    delta_squared = Fraction(3, 2)**dimension / expected
    rows.append(dict(dimension=dimension, normal_rank=rank,
                     primitive_normal_determinant=expected, saturation_index=index,
                     center_density_squared=str(delta_squared),
                     center_density=math.sqrt(float(delta_squared)),
                     status='Reproduces the reference bound; no improvement.',
                     minimum_dependency='Published minimum 6 of P48n.'))
(P/'normal_witness_checks.json').write_text(json.dumps(rows, indent=2))
print(json.dumps(rows, indent=2))
