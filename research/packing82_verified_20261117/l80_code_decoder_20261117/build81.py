from pathlib import Path
from fractions import Fraction as F
import json,math,numpy as np
from flint import fmpz_mat as Z,fmpq_mat as Q
p=Path(__file__).resolve().parent;a=np.load(p.parent/'bachoc_nebe80_reconstruction_20261117/lattice80.npz');G=Z(a['gram'].tolist());M=Z(a['twice_root_coordinate_basis'].tolist());E=Z(a['e8_twice_euclidean_basis'].tolist());H=np.concatenate([a['T1'],a['T2']]);case=json.loads((p/'case101.json').read_text());run=json.loads((p/'case101_contact.json').read_text());z=Q([case['diagonal_z_root'][run['glue']]]);row=[]
for state,y2 in zip(run['states'],run['local_points_times2']):
 s=Q([(np.array([(state>>i)&1 for i in range(8)])@H).tolist()]);y=Z([y2])*E.inv();assert all(x.denominator==1 for x in y.entries());row+=list((2*s+4*y+z).tolist()[0])
v=Q([row])*M.inv();assert all(x.denominator==1 for x in v.entries());r=Z([case['r']]);res=r-2*Z([[int(x)for x in v.entries()]]);assert (res*G*res.transpose())[0,0]==24;assert (r*G*r.transpose())[0,0]%4==0
# Use the short residual representative, differing from the frozen target by2L.
x=Q(res)/2;cross=Q(G)*x.transpose();Gram=Q(81,81)
for i in range(80):
 for j in range(80):Gram[i,j]=G[i,j]
 Gram[i,80]=Gram[80,i]=cross[i,0]
Gram[80,80]=8;assert Gram.det()==2
D=F(2)**80;assert F(2)**81/F(str(Gram.det()))==D

def at(q,N):return sum((F((-1)**j,(2*j+1)*q**(2*j+1))for j in range(N+1)),F(0))
pi=16*at(5,161)-4*at(239,160);vol2=F(2**82,math.prod(range(1,82,2))**2)*pi**80;floor=F('2.87777416803990e-17');assert floor**2<D*vol2
prior=F('79172296088172630531378644744173191017334418459007976066782137593751953125/69620512594693245636748678160327321048267519688704');external=F(7,2)**81/(14*41**12)
out={'dimension':81,'status':'PRIMARY_COMPLETE_CODE_EXCLUSION; independent audit pending','parent_minimum':8,'parent_determinant':1,'coset_depth':6,'vertical_squared_height':2,'gram_determinant':'2','center_density_squared':str(D),'ordinary_density_floor':str(floor),'short_coset_numerator':res.tolist()[0],'original_coset_numerator':r.tolist()[0],'contact_lattice_vector':[[str(t)for t in row]for row in v.tolist()][0],'gram':[[str(t)for t in row]for row in Gram.tolist()],'ratio_to_pinned_repo':math.sqrt(float(D/prior)),'ratio_to_derived_MW80_A1':math.sqrt(float(D/external)),'pinned_repo_center_density_squared':str(prior),'external_center_density_squared':str(external)};(p/'candidate81.json').write_text(json.dumps(out,indent=2,default=int));print({k:v for k,v in out.items()if k not in ['gram','short_coset_numerator','original_coset_numerator','contact_lattice_vector']})
