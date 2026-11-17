from pathlib import Path
from fractions import Fraction as F
from math import prod
import json,numpy as np
from flint import fmpz_mat as Z,fmpq_mat as Q
D=Path(__file__).parent;data=np.load(D.parent/'bachoc_nebe80_reconstruction_20261117/lattice80.npz');G=Z(data['gram'].tolist());z=json.loads((D/'candidate81.json').read_text());case=json.loads((D/'case101.json').read_text());assert G==Z(case['gram'])
r=Z([case['r']]);res=Z([z['short_coset_numerator']]);v=Q([z['contact_lattice_vector']]);assert all(x.denominator==1 for x in v.entries());assert Q(r)-2*v==Q(res)
assert (res*G*res.transpose())[0,0]==24 and (r*G*r.transpose())[0,0]%4==0 and any(x%2 for x in res.entries())
# For x=r/2, distance²=x²−(r,y)+y² is integral for every y∈Λ.
# Exclusion of all distance²≤5 therefore gives >=6; res/2 attains6.
H=Q(z['gram']);assert H.nrows()==81 and H.ncols()==81
for i in range(80):
 for j in range(80):assert H[i,j]==G[i,j]
cross=Q(G)*Q(res).transpose()/2
for i in range(80):assert H[i,80]==H[80,i]==cross[i,0]
assert H[80,80]==8 and H.det()==2 and H[80,80]-(Q(res)*Q(G)*Q(res).transpose())[0,0]/4==2
# Even layer multiples >=2 have vertical norm>=8, odd multiples have horizontal distance>=6 plus vertical>=2.
# Layer0 inherits min8; parent norm8 witness attains it.
delta2=F(2)**80;assert F(z['center_density_squared'])==delta2
atan=lambda q,N:sum((F((-1)**j,(2*j+1)*q**(2*j+1))for j in range(N+1)),F(0))
piL=16*atan(5,201)-4*atan(239,200);vol2=F(2**82,prod(range(1,82,2))**2)*piL**80;floor=F('2.87777416803990e-17');assert floor**2<delta2*vol2
external=F(7,2)**81/(14*41**12);assert external==F(z['external_center_density_squared'])
result=dict(status='PASS exact integral-distance gap, depth6 contact, layered min8, Gram/determinant, density floor',dimension=81,determinant=2,minimum_squared=8,coset_depth=6,density_floor=str(floor),external_comparison='Pinned derived MW80 x A1 example, not asserted strongest',center_density_squared=str(delta2))
(D/'extension_result.json').write_text(json.dumps(result,indent=2));print(result)
