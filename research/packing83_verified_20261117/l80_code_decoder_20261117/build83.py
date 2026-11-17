from pathlib import Path
from fractions import Fraction as F
import json,itertools,math,numpy as np
from flint import fmpz_mat as Z,fmpq_mat as Q,fmpq as q
p=Path(__file__).resolve().parent;parent=np.load(p.parent/'bachoc_nebe80_reconstruction_20261117/lattice80.npz');G=Z(parent['gram'].tolist());M=Z(parent['twice_root_coordinate_basis'].tolist());geo=json.loads((p/'geometry83.json').read_text());sources={ix:Z([json.loads((p/f'case{ix}.json').read_text())['r']])for ix in[101,108,203]};actions={}
def transfer(name,src):
 rec=geo[name];bp=rec['block_permutation'];assert sorted(bp)==list(range(10));P=Z(80,80)
 for i,j in enumerate(bp):
  for k in range(8):P[8*i+k,8*j+k]=1
 A=M*P*M.inv();assert all(x.denominator==1 for x in A.entries());A=Z([[int(x)for x in row]for row in A.tolist()]);assert abs(A.det())==1 and A*G*A.transpose()==G
 v=sources[src]*A;assert v.tolist()[0]==rec['lift'];r=Z([rec['class']]);assert all(int(x)%2==0 for x in(v-r).entries());actions[name]=A.tolist();return r
ra=transfer('a',101);rb=transfer('b',101);rc=transfer('c',108);rac=transfer('cross_ac',203);rbc=transfer('cross_bc',203)
for v,w in[(ra+rc,rac),(rb+rc,rbc)]:assert all(int(x)%2==0 for x in(v-w).entries())
for ix,v,cut,nodes in[(301,ra+rb,3,7416859),(302,ra+rb+rc,4,149419024)]:
 r=Z([json.loads((p/f'case{ix}.json').read_text())['r']]);assert all(int(x)%2==0 for x in(v-r).entries());assert(r*G*r.transpose())[0,0]%4==0;run=json.loads((p/f'case{ix}_run.json').read_text());assert run['status']=='EXHAUSTIVE_NO_VECTOR_PENDING_AUDIT' and run['nodes']==nodes and run['leaves']==0
V=Q([[2,0,q(3,4)],[0,2,q(3,4)],[q(3,4),q(3,4),2]]);assert V.det()==q(23,4);L=V-Q([[q(9,10)*int(i==j)for j in range(3)]for i in range(3)]);assert all(Q([row[:k]for row in L.tolist()[:k]]).det()>0 for k in[1,2,3]);mins=[None]*8
for v in itertools.product(range(-3,4),repeat=3):
 if not any(v):continue
 r=Q([list(v)]);val=(r*V*r.transpose())[0,0];par=sum((v[i]%2)<<i for i in range(3));mins[par]=val if mins[par]is None else min(mins[par],val)
assert mins==list(map(q,['8','2','2','4','2','5/2','5/2','3']))
# Outside the cube norm²>=16, so V>=9/10 I gives value>8.
R=Q([ra.tolist()[0],rb.tolist()[0],rc.tolist()[0]])/2;H=Q(83,83);C=R*G;T=R*G*R.transpose()+V
for i in range(80):
 for j in range(80):H[i,j]=G[i,j]
 for j in range(3):H[i,80+j]=H[80+j,i]=C[j,i]
for i in range(3):
 for j in range(3):H[80+i,80+j]=T[i,j]
assert H.det()==q(23,4);D=F(2)**85/23
at=lambda q,N:sum((F((-1)**j,(2*j+1)*q**(2*j+1))for j in range(N+1)),F(0));pi=16*at(5,161)-4*at(239,160);vol2=F(2**84,math.prod(range(1,84,2))**2)*pi**82
# Select a conservative15-significant-digit floor using float, then certify rationally.
fapprox=math.pi**41.5/math.gamma(42.5)*math.sqrt(float(D));power=math.floor(math.log10(fapprox))-14;floor=F(math.floor(fapprox/10**power))*F(10)**power
while floor**2>=D*vol2:floor-=F(10)**power
prior=F('154742504910672534362390528/125');external=F(7,2)**83/(1372*41**12)
out={'dimension':83,'status':'PRIMARY_COMPLETE_PARITY_EXTENSION; independent audit pending','minimum_squared':8,'gram_determinant':'23/4','vertical_gram':[[str(x)for x in row]for row in V.tolist()],'vertical_parity_minima':[str(x)for x in mins],'horizontal_parity_minima':['0','6','6','4','6','11/2','11/2','5'],'gram':[[str(x)for x in row]for row in H.tolist()],'coset_numerators':[ra.tolist()[0],rb.tolist()[0],rc.tolist()[0]],'automorphisms':actions,'center_density_squared':str(D),'ordinary_density_floor':str(floor),'density_float_approx':fapprox,'pinned_repo_center_density_squared':str(prior),'ratio_to_pinned_repo':math.sqrt(float(D/prior)),'external_MW80_A3_center_density_squared':str(external),'ratio_to_external_example':math.sqrt(float(D/external))};(p/'candidate83.json').write_text(json.dumps(out,indent=2,default=int));print({k:v for k,v in out.items()if k not in['gram','coset_numerators','automorphisms']})
