from pathlib import Path
from fractions import Fraction as F
import json
from flint import fmpq as Q,arb,arb_poly,ctx
if not __debug__:raise RuntimeError('Assertions required')
ctx.prec=384;D=Path(__file__).resolve().parent;s=json.loads((D/'witness.json').read_text());rad=json.loads((D/'radial_exact.json').read_text());b=next(x for x in json.loads((D/'budget.json').read_text())['dimensions']if x['dimension']==5)
assert b['geometry_mode']=='spherical_transverse' and F(b['short_radius'])==F(27,25)
assert list(map(F,b['physical_function']['exact_rational_coefficients']))==[F(str(x))for x in s['x']]
assert b['physical_function']['additional_physical_terms']=='none'
a=lambda v:arb(Q(str(v)))
P=arb_poly([a(F(x))for x in rad['P']]);P1=P.derivative()-P;P2=P1.derivative()-P1
bins=[(F(x['lo']),F(x['hi']),int(x['cap']),F(x['eta']))for x in b['bins']];assert len(bins)==132
assert bins[0][0]==F(3,2) and bins[-1][1]==F(54,25)
assert all(x[1]==y[0]for x,y in zip(bins,bins[1:]))
for j,bd in enumerate(b['boundaries']):
 assert F(bd['base'])==(bins[j][0]if j<len(bins)else bins[-1][1])
 assert F(bd['eta'])==min(bins[j-1][3]if j else F(0),bins[j][3]if j<len(bins)else F(0))
caps=json.loads((D/'cap_audit.json').read_text());assert caps['pass'];table=[(F(z['lower_base']),z['cap'])for z in caps['table']]
pi=arb.pi();leaves=0;maxdepth=0
for lo,hi,cap,eta in bins:
 rho=F(27,25)**2-lo*lo/4
 assert lo>=F(3,2) and lo<hi<=F(54,25) and eta>=0 and cap>=1
 assert cap==[v for c,v in table if c<=lo][-1]
 # Verify on the larger X interval obtained from interval pi, using derivatives
 # with respect to X rather than the source verifier's radius derivatives.
 xlo=pi*a(lo*lo);xhi=pi*a(hi*hi);I=xlo.union(xhi);stack=[(I,0)]
 while stack:
  I,depth=stack.pop();m=I.mid();h=I.rad();g=P(m)*(-m).exp()+a(cap*eta);gp=P1(m)*(-m).exp()
  bound2=sum((abs(c)for c in P2(arb_poly([m,h])).coeffs()),arb(0))*(-I.lower()).exp()
  assert g+abs(gp)*h+bound2*h*h/2<0 or depth<30
  if g+abs(gp)*h+bound2*h*h/2<0:leaves+=1;maxdepth=max(maxdepth,depth);continue
  stack.extend([(I.lower().union(m),depth+1),(m.union(I.upper()),depth+1)])
out={'pass':True,'closed_bins':len(bins),'leaves':leaves,'maxdepth':maxdepth,'caps':'independently verified spherical, Rankin and centroid bounds','boundary_rule':'minimum adjacent allowances, zero at support endpoints'}
(D/'budget_audit.json').write_text(json.dumps(out,indent=2));print(out)
