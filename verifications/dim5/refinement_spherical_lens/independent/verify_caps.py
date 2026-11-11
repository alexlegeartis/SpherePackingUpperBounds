"""Independent exact 5D transverse lens caps, including explicit S^3 polynomials."""
from pathlib import Path
from fractions import Fraction as F
import json,math
from flint import fmpq as Q,fmpq_poly as P
if not __debug__:raise RuntimeError('Assertions required')
D=Path(__file__).resolve().parent;b=next(x for x in json.loads((D/'budget.json').read_text())['dimensions']if x['dimension']==5)
assert b['geometry_mode']=='spherical_transverse';r=F(27,25);A=(r*r-1)**2;cmax=F(41,25);X=P([0,1]);table=[]
rows=b['exact_cap_geometry'];assert list(map(lambda z:F(z['lower_base']),rows))==list(map(F,['3/2','31/20','8/5','41/25','17/10','7/4','9/5','46/25','48/25']))
for i,row in enumerate(rows):
 c=F(row['lower_base']);cap=int(row['cap']);U=r*r-c*c/4;delta=A/(c*c)
 if i<3:
  s=row['spherical_certificate'];L=1-c*c/4-delta/4;angle=1-(1-delta)/(2*U)
  assert s['packing_dimension']==5 and s['ambient_dimension']==4 and s['projection']=='transverse' and F(s['base_lower'])==c
  assert F(s['norm_squared_lower'])==L>0 and F(s['norm_squared_upper'])==U>0
  assert F(s['transverse_longitudinal_difference_squared'])==delta
  assert F(s['max_inner_product'])==angle<1 and F(s['monotonicity_margin'])==L-U+1-delta>0
  assert F(s['valid_through_base'])==cmax
  lower=1-cmax*cmax/4-A/(4*cmax*cmax);assert F(s['uniform_norm_squared_lower'])==lower>0
  assert c**4>A and F(s['angular_monotonicity_positive_remainder'])==A*(4*r*r-A)>0
  cs=list(map(Q,s['gegenbauer_coefficients']));assert cs[0]==1 and min(cs)>=0
  # Dimension4 normalized Gegenbauer is Chebyshev U_l/(l+1).
  pol=P([])
  for k,co in enumerate(cs):pol+=co*sum(((-1)**j*math.comb(k-j,j)*(2*X)**(k-2*j)for j in range(k//2+1)),P([]))/Q(k+1)
  assert pol==P(list(map(Q,s['polynomial']))) and pol(1)==sum(cs)==Q(s['rational_bound'])
  assert cap==s['integer_bound'] and cap<=pol(1)<cap+1
  lo=Q(-1);hi=Q(str(angle));assert pol(lo)<0 and pol(hi)<0;seq=[pol,pol.derivative()]
  while seq[-1].degree()>0:
   z=-(seq[-2]%seq[-1]);assert z;seq.append(z)
  def V(t):
   ss=[int(p(t)>0)-int(p(t)<0)for p in seq if p(t)];return sum(a!=b for a,b in zip(ss,ss[1:]))
  assert V(lo)==V(hi)
 else:
  assert 0<=U<F(1,2);centroid=int(1/(1-2*U));dot=U-(1-delta)/2
  assert dot<0 and F(row['obtuse_vector_dot_upper'])==dot
  assert row['obtuse_vector_dimension']==4 and row['obtuse_cap']==5
  assert row['centroid_cap']==centroid and cap==min(centroid,5)
 table.append({'lower_base':str(c),'cap':cap})
assert [z['cap']for z in table]==[10,9,8,5,5,5,3,2,1]
(D/'cap_audit.json').write_text(json.dumps({'pass':True,'table':table,'scope':'Explicit spherical polynomials/global Sturm signs, rational projection/monotonicity margins, elementary Rankin and centroid caps'},indent=2));print('PASS exact transverse spherical10/9/8, Rankin5 and centroid caps')
