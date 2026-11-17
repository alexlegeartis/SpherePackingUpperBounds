from pathlib import Path
import json,hashlib,struct,sys,math
from fractions import Fraction as F
import numpy as np
from flint import fmpq_mat as M,fmpq as Q,fmpz_mat as Z
if not __debug__:raise RuntimeError('Assertions required')
P=Path(__file__).resolve().parent/'input'
for f,h in json.load(open(P/'SHA256SUMS.json')).items():assert hashlib.sha256((P/f).read_bytes()).hexdigest()==h,f
E=json.load(open(P/'random56_37_exact.json'));G=M(E['G']);r=list(map(int,E['target_numerator']));L=M(E['mu_rational']);ds=list(map(Q,E['diagonal_rational']));n=56
assert len(r)==n and E['denominator']==3 and E['threshold_numerator']==54
assert all(0<=x<=2 for x in r) and any(r)
assert all(ds[i]>0 and L[i,i]==1 and all(L[i,j]==0 for j in range(i+1,n)) for i in range(n))
D=M([[ds[i] if i==j else 0 for j in range(n)]for i in range(n)])
assert L*D*L.transpose()==G==G.transpose()
assert all(G[i,j].denominator==1 for i in range(n)for j in range(n));assert all(int(G[i,i])%2==0 for i in range(n))
assert sum(r[i]*int(G[i,j])*r[j] for i in range(n)for j in range(n))%6==0
# Independently parse every binary64 endpoint by bit pattern and compare exact fractions.
tokens=list(map(int,(P/'random56_37_bits.txt').read_text().split()));assert tokens[:2]==[56,54];pos=2
assert tokens[pos:pos+n]==r;pos+=n
assert tokens[pos:pos+n*n]==[int(G[i,j])for i in range(n)for j in range(n)];pos+=n*n
for val in ds+[L[i,j]for i in range(n)for j in range(n)]:
 lo,hi=[struct.unpack('d',struct.pack('Q',tokens[pos+k]))[0]for k in range(2)];pos+=2
 assert math.isfinite(lo) and math.isfinite(hi) and F.from_float(lo)<=F(str(val))<=F.from_float(hi)
assert pos==len(tokens)
B=json.load(open(P/'cvp_basis_actions.json'));T=Z(B['T']);base=np.load(P/'parent/parent56.npz');H=Z(base['G'].tolist());assert abs(T.det())==1 and T*H*T.transpose()==G
C=json.load(open(P/'random37_contact.json'));v=C['residual_numerator'];assert C['denominator']==3 and all((v[i]-r[i])%3==0 for i in range(n));assert sum(v[i]*int(G[i,j])*v[j] for i in range(n)for j in range(n))==60
parentv=list((Z([v])*T).tolist()[0]);assert parentv==C['parent_coordinate_numerator']
assert H.det()==29**8
# Verify literal block Grams and exact Schur determinant.
x=M([[Q(t,3)for t in parentv]]);HH=M(H);cross=HH*x.transpose();xx=(x*HH*x.transpose())[0,0];assert xx==Q(20,3)
for rec in json.load(open(P/'extension_candidates.json')):
 dim=rec['dimension'];A=M(rec['gram']);assert dim in [57,58] and A.nrows()==dim and A==A.transpose()
 for i in range(56):
  for j in range(56):assert A[i,j]==HH[i,j]
  assert A[i,56]==cross[i,0]
  if dim==58:assert A[i,57]==-cross[i,0]
 assert A[56,56]==10
 if dim==58:assert A[57,57]==10 and A[56,57]==Q(-5)
 det=Q(10,3)*29**8 if dim==57 else Q(25,3)*29**8
 assert A.det()==det==Q(rec['gram_determinant'])
 cd=Q(5,2)**dim/det;assert cd==Q(rec['center_density_squared'])
 print(dim,'EXACT GRAM/DETERMINANT/CENTER DENSITY PASS',str(det),flush=True)
print('HASHES, LDL, BINARY ENCLOSURES, UNIMODULAR BASIS, CONGRUENCE AND CONTACT PASS')
