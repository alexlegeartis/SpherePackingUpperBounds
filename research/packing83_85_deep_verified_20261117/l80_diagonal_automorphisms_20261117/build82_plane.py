from pathlib import Path
from fractions import Fraction as F
import json,math
from flint import fmpq_mat as Q
import numpy as np
p=Path(__file__).resolve().parent;parent=np.load(p.parent/'bachoc_nebe80_reconstruction_20261117/lattice80.npz');G=Q(parent['gram'].tolist());items=json.loads((p/'deep_plane.json').read_text());R=Q([x['class']for x in items[:2]])/2;V=Q([[2,-1],[-1,2]]);H=Q(82,82);C=R*G;T=R*G*R.transpose()+V
for i in range(80):
 for j in range(80):H[i,j]=G[i,j]
 for j in range(2):H[i,80+j]=H[80+j,i]=C[j,i]
for i in range(2):
 for j in range(2):H[80+i,80+j]=T[i,j]
assert H.det()==3
# A2 minimum2; twice every nonzero vector has norm>=8. Every nonzero parity class has horizontal depth>=6.
D=F(2)**82/3
at=lambda q,N:sum((F((-1)**j,(2*j+1)*q**(2*j+1))for j in range(N+1)),F(0));pi=16*at(5,201)-4*at(239,200);vol2=pi**82/F(math.factorial(41))**2
fapprox=math.pi**41/math.factorial(41)*math.sqrt(float(D));power=math.floor(math.log10(fapprox))-14;floor=F(math.floor(fapprox/10**power))*F(10)**power
while floor**2>=D*vol2:floor-=F(10)**power
prior=F(2)**86/55;ext=F(7,2)**82/(147*41**12)
z={'dimension':82,'minimum_squared':8,'gram_determinant':'3','vertical_gram':[['2','-1'],['-1','2']],'center_density_squared':str(D),'ordinary_density_floor':str(floor),'gram':[[str(x)for x in row]for row in H.tolist()],'ratio_to_current':math.sqrt(float(D/prior)),'external_MW80_A2_center_density_squared':str(ext),'ratio_to_external_example':math.sqrt(float(D/ext))}
(p/'candidate82_plane.json').write_text(json.dumps(z,indent=2));print({k:v for k,v in z.items()if k!='gram'})
