"""Reconstruct the frozen degree-16 polynomial from its numerical checkpoint."""
from pathlib import Path
from decimal import Decimal,ROUND_HALF_EVEN,ROUND_CEILING
from fractions import Fraction as Q
import json

HERE=Path(__file__).resolve().parent
if not __debug__:raise RuntimeError('assertions must be enabled')
PACKAGE=HERE.parents[1]/'verifications/dim4/refinement_energy20261115'
checkpoint=json.loads((HERE/'checkpoint.json').read_text(),parse_float=str)
reconstructed=[]
for d,value in enumerate(checkpoint['coefficients']):
    x=Decimal(value)
    if x<0 and d>=11:q=Q(0)
    else:
        rounded=(x*10**9).to_integral_value(rounding=ROUND_CEILING if x<0 else ROUND_HALF_EVEN)
        q=Q(int(rounded),10**9)
    reconstructed.append(q*Q(999,1000)*Q(1,10**6))
published=json.loads((PACKAGE/'polynomials.json').read_text())
assert reconstructed==list(map(Q,published['1']['coefficients']))
print('PASS: exact frozen-checkpoint export matches every published coefficient')
