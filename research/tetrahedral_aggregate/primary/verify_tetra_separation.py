"""Exact two-point obstruction to any affine envelope for the verified 22D Q."""
import json
from pathlib import Path
from fractions import Fraction as F
root=Path(__file__).resolve().parent
r=json.loads((root/'tetra_full_geometry22_check.json').read_text())
assert r['dimension']==22 and F(r['q'])==F(21,20)
q=F(r['q']);lo=1-q/(2*(1-q/4));hi=1-1/(2*(q-F(1,4)))
assert lo==F(17,59) and hi==F(3,8)
cs=list(map(F,r['polynomial']))
def value(t):
 p=F(0)
 for c in reversed(cs):p=p*t+c
 return p
lower=value(hi)/(hi-F(1,3));upper=value(lo)/(lo-F(1,3))
assert lower==F(7034984040475208405061,65536000000000000000)
assert upper==F(1191611909750310730113443859,15728334578000000000000000)
assert lower>upper
print('No affine slope exists: a >=',lower,'but a <=',upper)
print('PASS: strict feasible-cone separation, conditional on the separately verified tetrahedron signs.')
