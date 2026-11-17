from pathlib import Path
from fractions import Fraction as F
import json,numpy as np,math
from flint import fmpz_mat as Z,fmpq_mat as Q,fmpq as q
p=Path(__file__).resolve().parent;a=np.load(p.parent/'bachoc_nebe80_reconstruction_20261117/lattice80.npz');G=Z(a['gram'].tolist());M=Z(a['twice_root_coordinate_basis'].tolist());c101=json.loads((p/'case101.json').read_text());c108=json.loads((p/'case108.json').read_text());c203=json.loads((p/'case203.json').read_text());b81=json.loads((p/'candidate81.json').read_text());rec=c203['second_coset'];bp=rec['block_permutation'];P=Z(80,80)
assert sorted(bp)==list(range(10))
for i,j in enumerate(bp):
 for k in range(8):P[8*i+k,8*j+k]=1
A=M*P*M.inv();assert all(x.denominator==1 for x in A.entries());A=Z([[int(x)for x in row]for row in A.tolist()]);assert abs(A.det())==1 and A*G*A.transpose()==G
s=Z([c108['r']])*A;assert s.tolist()[0]==rec['lift'];r2=Z([rec['class']]);assert all(int(x)%2==0 for x in(s-r2).entries());r1=Z([c101['r']]);r3=Z([c203['r']]);assert all(int(x)%2==0 for x in(r1+r2-r3).entries());assert(r3*G*r3.transpose())[0,0]%4==2
assert(r1*G*r1.transpose())[0,0]%4==0 and(s*G*s.transpose())[0,0]%4==0
for ix,nodes,leaves in [(108,5099300880,100),(203,2489483280,0)]:
 run=json.loads((p/f'case{ix}_run.json').read_text());assert run['status']=='EXHAUSTIVE_NO_VECTOR_PENDING_AUDIT' and run['nodes']==nodes
# Height forms: parity10/01>=2, parity11>=5/2, nonzero even>=8.
V=Q([['2','3/4'],['3/4','2']]);assert V.det()==q('55/16');assert V[0,0]-q('5/4')>0 and (V-Q([[q('5/4'),0],[0,q('5/4')]])).det()==0
R=Q([b81['short_coset_numerator'],rec['class']])/2;H=Q(82,82);C=R*G;T=R*G*R.transpose()+V
for i in range(80):
 for j in range(80):H[i,j]=G[i,j]
 for j in range(2):H[i,80+j]=H[80+j,i]=C[j,i]
for i in range(2):
 for j in range(2):H[80+i,80+j]=T[i,j]
assert H.det()==q('55/16');D=F(2)**86/55
# Rational decimal lower floor selected conservatively, independently checked below.
density='8.5669123663438021942076969170632968764e-18';floor=F('8.56691236634380e-18')
at=lambda q,N:sum((F((-1)**j,(2*j+1)*q**(2*j+1))for j in range(N+1)),F(0));pi=16*at(5,161)-4*at(239,160);assert floor**2<D*pi**82/math.factorial(41)**2
prior=F('6185335631888486760263956620638530548229251442109998130217354499511871337890625/5714934664858416272284567663612103329623072090271776768');external=F(7,2)**82/(147*41**12)
out={'dimension':82,'status':'PRIMARY_COMPLETE_MIXED_DEPTH_CONSTRUCTION; independent audit pending','parent_minimum':8,'minimum_squared':8,'half_coset_depth_lower_bounds':['6','6','11/2'],'vertical_gram':[['2','3/4'],['3/4','2']],'gram_determinant':'55/16','gram':[[str(x)for x in row]for row in H.tolist()],'coset_numerators':[b81['short_coset_numerator'],rec['class']],'coordinate_automorphism':A.tolist(),'center_density_squared':str(D),'ordinary_density_floor':str(floor),'ordinary_density_decimal_approx':str(density),'pinned_repo_center_density_squared':str(prior),'ratio_to_pinned_repo':math.sqrt(float(D/prior)),'external_MW80_A2_center_density_squared':str(external),'ratio_to_external_example':math.sqrt(float(D/external))};(p/'candidate82.json').write_text(json.dumps(out,indent=2,default=int));print({k:v for k,v in out.items()if k not in ['gram','coordinate_automorphism','coset_numerators']})
