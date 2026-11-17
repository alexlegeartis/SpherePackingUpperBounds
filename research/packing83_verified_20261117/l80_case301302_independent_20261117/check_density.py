from pathlib import Path
from fractions import Fraction as F
from math import prod
import json,decimal
from flint import fmpq_mat as Q
D=Path(__file__).parent;z=json.loads((D/'candidate83.json').read_text());H=Q(z['gram']);ind=Q(json.loads((D.parent/'l80_83_geometry_independent_20261117/gram83.json').read_text()));assert H==ind and str(H.det())=='23/4'
delta2=F(2)**85/23;assert delta2==F(z['center_density_squared']);floor=F('2.56962386473426e-18');assert floor==F(z['ordinary_density_floor'])
at=lambda q,N:sum((F((-1)**j,(2*j+1)*q**(2*j+1))for j in range(N+1)),F(0));piL=16*at(5,201)-4*at(239,200);assert floor**2<delta2*F(2**84,prod(range(1,84,2))**2)*piL**82
external=F(7,2)**83/(4*7**3*41**12);assert external==F(z['external_MW80_A3_center_density_squared']);decimal.getcontext().prec=60;gain=(decimal.Decimal((delta2/external).numerator)/decimal.Decimal((delta2/external).denominator)).sqrt()*100-100
out={'status':'PASS matching independently reconstructed full Gram and rigorous density floor','density_floor':'2.56962386473426e-18','external_gain_percent':str(gain),'external_comparator':'Pinned derived MW80 x rescaled A3, not strongest-prior claim'};(D/'density_result.json').write_text(json.dumps(out,indent=2));print(out)
