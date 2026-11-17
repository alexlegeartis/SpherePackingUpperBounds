from pathlib import Path
from fractions import Fraction as F
from math import factorial
import json,numpy as np
from flint import fmpq_mat as Q
D=Path(__file__).resolve().parent;z=json.loads((D/'candidate82.json').read_text());G=Q(np.load(D.parent/'bachoc_nebe80_reconstruction_20261117/lattice80.npz')['gram'].tolist());H=Q(z['gram']);c=json.loads((D/'case203.json').read_text());old=json.loads((D.parent/'l80_case101_independent_20261117/case101.json').read_text());R=Q(z['coset_numerators'])/2
for row,bits in zip(z['coset_numerators'],[old['r'],c['second_coset']['class']]):assert all((int(x)-y)%2==0 for x,y in zip(row,bits))
assert H.nrows()==H.ncols()==82
for i in range(80):
 for j in range(80):assert H[i,j]==G[i,j]
 for j in range(2):assert H[i,80+j]==H[80+j,i]==(G*R.transpose())[i,j]
V=Q([[H[80+i,80+j]for j in range(2)]for i in range(2)])-R*G*R.transpose();assert V==Q([['2','3/4'],['3/4','2']]);assert H.det()==V.det() and str(V.det())=='55/16'
# Q >=5/4(m²+n²). If m²+n²>=7 then Q>8. Exhaust all smaller vectors, including each parity minimum.
mins={}
for m in range(-3,4):
 for n in range(-3,4):
  if (m,n)==(0,0):continue
  q=F(2)*(m*m+n*n)+F(3,2)*m*n;assert q>=F(5,4)*(m*m+n*n)
  k=(m%2,n%2);mins[k]=min(mins.get(k,q),q)
assert mins=={(1,1):F(5,2),(1,0):F(2),(0,1):F(2),(0,0):F(8)}
for k,depth in [((1,0),6),((0,1),6),((1,1),F(11,2)),((0,0),0)]:assert mins[k]+depth>=8
# Origin layer uses parent minimum8. Parent norm8 witness attains the final minimum.
delta2=F(2)**82/F(55,16);assert delta2==F(2)**86/55==F(z['center_density_squared'])
atan=lambda q,N:sum((F((-1)**j,(2*j+1)*q**(2*j+1))for j in range(N+1)),F(0))
piL=16*atan(5,201)-4*atan(239,200);floor=F('8.56691236634380e-18');assert floor**2<delta2*piL**82/factorial(41)**2
# MW80 min14 det41^12 times A2 rescaled to min14 has determinant147*41^12.
external=F(7,2)**82/(147*41**12);assert external==F(z['external_MW80_A2_center_density_squared'])
import decimal
decimal.getcontext().prec=60;dec=lambda x:decimal.Decimal(x.numerator)/decimal.Decimal(x.denominator);gain=100*((dec(delta2/external)).sqrt()-1)
out={'status':'PASS exact Gram, parity minima and density conditional on completed frozen exclusions','dimension':82,'minimum_squared':8,'determinant':'55/16','density_floor':str(floor),'center_density_squared':str(delta2),'external_comparator':'Pinned derived MW80 x A2, not strongest-previous claim','external_gain_percent':str(gain)};(D/'extension_result.json').write_text(json.dumps(out,indent=2));print(out)
