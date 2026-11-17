from pathlib import Path
from fractions import Fraction as F
import json,struct,math,hashlib
from flint import fmpq_mat as M,fmpq as Q,fmpz_mat as Z
P=Path(__file__).resolve().parent; old=P.parent/'cm57_independent_20261117/input'; src=P.parent/'cm57_deeper_20261117'
assert hashlib.sha256((P/'input.txt').read_bytes()).hexdigest()=='ad7610e62d751424f833dd8ab779f3153363c10a55c413eeef29454bd9ee7cf1'
E=json.loads((old/'random56_37_exact.json').read_text());G=M(E['G']);L=M(E['mu_rational']);ds=list(map(Q,E['diagonal_rational']));n=56
assert L*M([[ds[i] if i==j else 0 for j in range(n)]for i in range(n)])*L.transpose()==G
assert all(ds[i]>0 and L[i,i]==1 and all(L[i,j]==0 for j in range(i+1,n)) for i in range(n))
t=list(map(int,(P/'input.txt').read_text().split()));assert t[:2]==[56,24];r=t[2:58];assert set(r)<={0,1} and any(r);pos=58
assert t[pos:pos+n*n]==[int(G[i,j])for i in range(n)for j in range(n)];pos+=n*n
for val in ds+[L[i,j]for i in range(n)for j in range(n)]:
 lo,hi=[struct.unpack('d',struct.pack('Q',t[pos+k]))[0]for k in range(2)];pos+=2
 assert math.isfinite(lo) and math.isfinite(hi) and F.from_float(lo)<=F(str(val))<=F.from_float(hi)
assert pos==len(t)
assert all(G[i,j].denominator==1 for i in range(n)for j in range(n))
R=M([r]);assert int((R*G*R.transpose())[0,0])%4==0
A=json.loads((old/'cvp_basis_actions.json').read_text());T=Z(A['T']);S=Z(A['old_to_new']);base=json.loads((old/'parent/primal_cm.json').read_text());C=S.inv()*Z(base['CM'])*S
I=Z([[int(i==j)for j in range(n)]for i in range(n)])
assert all(v.denominator==1 for v in C.entries()) and M(I)+C+C*C==M(n,n) and C*G*C.transpose()==G
assert abs(T.det())==1
contact=json.loads((src/'order2_527_contact.json').read_text());v=M([contact['residual_numerator']]);assert all((int(v[0,i])-r[i])%2==0 for i in range(n));assert (v*G*v.transpose())[0,0]==28
assert (v*C*G*v.transpose())[0,0]==-14
H=T.inv()*G*T.transpose().inv();assert H.det()==29**8
x=v*T/2;y=v*C*T/2
assert list((v*T).tolist()[0])==contact['parent_coordinate_numerator']
def at(q,N):return sum((F((-1)**k,(2*k+1)*q**(2*k+1))for k in range(N+1)),F(0))
pi=16*at(5,201)-4*at(239,200)
for rec in json.loads((src/'binary_extension_candidates.json').read_text()):
 d=rec['dimension'];B=M(rec['gram']);assert B==B.transpose()
 for i in range(n):
  for j in range(n):assert B[i,j]==H[i,j]
  assert B[i,56]==(x*H)[0,i]
  if d==58:assert B[i,57]==(y*H)[0,i]
 assert B[56,56]==10
 if d==58:assert B[57,57]==10 and B[56,57]==-2
 det=Q(3)*29**8 if d==57 else Q(27,4)*29**8
 assert B.det()==det==Q(rec['determinant'])
 delta2=F(5,2)**d/F(str(det));assert delta2==F(rec['center_density_squared'])
 vol2=pi**d/F(math.factorial(d//2)**2) if d%2==0 else F(2**((d+1)//2),math.prod(range(1,d+1,2)))**2*pi**(d-1)
 floor=F('1.61599884109140e-11' if d==57 else '5.58243739609337e-12');assert floor**2<delta2*vol2
print('PASS exact LDL/enclosures, CM, mod4, contact28, literal extension Grams/determinants, independent Machin density floors. Exclusion still separate.')
