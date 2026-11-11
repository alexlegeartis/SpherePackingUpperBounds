from fractions import Fraction as F
from pathlib import Path
import json, math
D=Path(__file__).parent
S=D/'continuum_three_point'
r=json.loads((S/'adaptive_cyclic40_radial111.json').read_text()); q=json.loads((S/'cyclic_angular40_exact.json').read_text())
p=list(map(F,r['physical_polynomial'])); B=F(q['B_exact']); N48=F(r['fourier_polynomial'][0]); M=F(r['neighbor_cost'])/B*F(1000001,1000000); pc=list(map(F,q['polynomial_coefficients']))
bench={25:.0013841907222857,26:.0009910238892216,27:.0007082297958617,28:.0005052542161057,29:.0003598581852089,30:.0002559028743732,31:.0001817083813917,32:.0001288432887595,33:.0000912356039023,34:.0000645221967438,35:.0000455743843107,36:.0000321530553313,37:.0000226586900106,38:.0000159506499105,39:.0000112168687009,40:.0000078801051697,41:.0000055306464395,42:.0000038780970907,43:.0000027169074727,44:.0000019017703144,45:.0000013300905665,46:.0000009295151556,47:.0000006490757338,48:.0000004529067791}
res=[]
for n in range(25,49):
 rising=F(1);N=F(0)
 for k,x in enumerate(p):
  N+=x*rising;rising*=F(n,2)+k
 mom=F(1);b0=F(0)
 for j in range((len(pc)+1)//2):
  b0+=pc[2*j]*mom;mom*=F(2*j+1,n+2*j)
 obj=(p[0]+M*B/b0)/N
 density=float(obj)*math.pi**(n/2)*1.5**(n/2)/math.gamma(n/2+1)
 out={'n':n,'N_ratio':float(N/N48),'b0':float(b0),'B_n':float(B/b0),'objective':float(obj),'density':density,'reference':bench[n],'ratio':density/bench[n]}
 print(out,flush=True);res.append(out)
(D/'restriction_screen.json').write_text(json.dumps(res,indent=2))
