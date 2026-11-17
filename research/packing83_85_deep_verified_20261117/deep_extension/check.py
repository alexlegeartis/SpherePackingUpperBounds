from pathlib import Path
import json,itertools,math,numpy as np
from fractions import Fraction as F
from flint import fmpq_mat as Q
p=Path(__file__).resolve().parent;root=p.parent;z=np.load(root/'bachoc_nebe80_reconstruction_20261117/lattice80.npz');G=Q(z['gram'].tolist());M=Q(z['twice_root_coordinate_basis'].tolist());Mi=M.inv();iso=json.load(open(root/'l80_diagonal_automorphisms_20261117/isometries.json'))['isometries'];results=[]
for n in [83,84,85]:
 c=json.load(open(p/f'candidate{n}.json'));rs=[x['class']for x in c['generators']];k=n-80;depth=[0]*2**k
 for mask in range(1,2**k):
  r=[sum(rs[j][i]for j in range(k)if mask>>j&1)%2 for i in range(80)];assert any(r)
  if str(mask) in c['known']:
   rec=c['known'][str(mask)];U=iso[rec['isometry_index']];bp=rec['block_permutation'];assert sorted(bp)==list(range(10));T=Q(80,80)
   for b in range(10):
    for i in range(8):
     for j in range(8):T[8*b+i,8*bp[b]+j]=U[i][j]
   A=M*T*Mi;assert all(x.denominator==1 for x in A.entries()) and A*G*A.transpose()==G and abs(A.det())==1
   src=json.load(open(root/f'l80_code_decoder_20261117/case{rec["source_case"]}.json'))['r'];lift=[int(x)for x in (Q([src])*A).entries()];assert lift==rec['lift'] and [x%2 for x in lift]==r;depth[mask]=6
  else:
   assert int((Q([r])*G*Q([r]).transpose())[0,0])%4==0
   old=list(map(int,(root/'l80_code_decoder_20261117/case101.txt').read_text().split()));t=list(map(int,(p/f'input{n}_{mask}.txt').read_text().split()));st=1+9*t[0]+64+256;X=(np.array(r,dtype=np.int64)@z['twice_root_coordinate_basis']).reshape(10,8)@z['e8_twice_euclidean_basis'];old[st:st+80]=X.ravel().tolist();assert t==old
   result=json.load(open(p/f'result{n}_{mask}.json'));assert result['status']=='EXHAUSTIVE_NO_VECTOR_PENDING_AUDIT';depth[mask]=4
 V=Q(c.get('vertical_gram',[[1+int(i==j)for j in range(k)]for i in range(k)]));det=int(V.det());assert det>0;mins=[10**6]*2**k
 for v in itertools.product(range(-5,6),repeat=k):
  if not any(v):continue
  val=int((Q([list(v)])*V*Q([list(v)]).transpose())[0,0]);m=sum((v[i]%2)<<i for i in range(k));mins[m]=min(mins[m],val)
 assert all(mins[m]+depth[m]>=8 for m in range(2**k))
 assert all(8*V.inv()[i,i] < 36 for i in range(k)) # q<8 forces |z_i|<6
 R=Q(rs)/2;H=Q(n,n);C=R*G;B=R*G*R.transpose()+V
 for i in range(80):
  for j in range(80):H[i,j]=G[i,j]
  for j in range(k):H[i,80+j]=H[80+j,i]=C[j,i]
 for i in range(k):
  for j in range(k):H[80+i,80+j]=B[i,j]
 assert H.det()==det
 D=F(2)**n/det;at=lambda q,N:sum((F((-1)**j,(2*j+1)*q**(2*j+1))for j in range(N+1)),F(0));pi=16*at(5,161)-4*at(239,160)
 if n%2:vol=F(2)**(n+1)*math.factorial((n+1)//2)*pi**((n-1)//2)/math.factorial(n+1)
 else:vol=pi**(n//2)/math.factorial(n//2)
 dd=D*vol**2;approx=math.sqrt(float(dd));power=math.floor(math.log10(approx))-14;step=F(10)**power;floor=math.floor(approx/float(step))*step
 while floor*floor>=dd:floor-=step
 out={'dimension':n,'determinant':det,'minimum_squared_at_least':8,'density_floor':str(floor),'density_approx':approx,'center_density_squared':str(D),'parity_minima':mins,'depths':depth,'gram':[[str(x)for x in row]for row in H.tolist()]};(p/f'verified{n}.json').write_text(json.dumps(out,indent=2));print(n,det,approx,'floor',floor,flush=True)
