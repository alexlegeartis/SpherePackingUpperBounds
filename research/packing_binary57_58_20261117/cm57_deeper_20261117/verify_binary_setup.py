from pathlib import Path
from fractions import Fraction as F
import json,math,struct
from flint import fmpz_mat,fmpq_mat
p=Path(__file__).parent;src=p.parent/'cm57_release_20261117';a=json.load(open(src/'cvp_basis_actions.json'));G=fmpz_mat(a['G']);C=fmpz_mat(a['old_to_new']);old=json.load(open(src/'parent/primal_cm.json'));CM=C.inv()*fmpz_mat(old['CM'])*C;assert all(x.denominator==1 for x in CM.entries());CM=fmpz_mat([[int(x)for x in row]for row in CM.tolist()]);I=fmpz_mat([[int(i==j)for j in range(56)]for i in range(56)]);assert I+CM+CM*CM==fmpz_mat(56,56) and CM*G*CM.transpose()==G
lines=(p/'order2_527_bits.txt').read_text().splitlines();original=(src/'random56_37_bits.txt').read_text().splitlines();assert lines[0]=='56' and lines[1]=='24' and lines[3:]==original[3:]
r=list(map(int,lines[2].split()));assert all(x in [0,1]for x in r) and any(r);R=fmpz_mat([r]);assert (R*G*R.transpose())[0,0]%4==0
contact=json.load(open(p/'order2_527_contact.json'));v=fmpz_mat([contact['residual_numerator']]);assert all((v[0,i]-r[i])%2==0 for i in range(56)) and (v*G*v.transpose())[0,0]==28
s=[int(x)%2 for x in (R*CM).entries()];assert s!=r and any(s)
run=json.load(open(p/'order2_527_complete.json'));assert run['status']=='EXHAUSTIVE_NO_VECTOR' and run['nodes']<10000000000
# Rational Machin lower enclosure, same elementary alternating-series identity.
def at(q,N):return sum((F((-1)**k,(2*k+1)*q**(2*k+1))for k in range(N+1)),F(0))
pi=16*at(5,121)-4*at(239,120);result=[]
for a in json.load(open(p/'binary_extension_candidates.json')):
 n=a['dimension'];det=F(str(fmpq_mat(a['gram']).det()));assert det==F(a['determinant']);delta2=F(5,2)**n/det;assert delta2==F(a['center_density_squared'])
 if n%2==0:vol2=pi**n/F(math.factorial(n//2)**2)
 else:vol2=F(2**((n+1)//2),math.prod(range(1,n+1,2)))**2*pi**(n-1)
 floor=F('1.61599884109140e-11' if n==57 else '5.58243739609337e-12');assert floor**2<delta2*vol2
 result.append({'dimension':n,'rigorous_density_floor':str(floor),'status':'Primary complete exclusion and exact construction pass; independent audit pending'})
(p/'binary_setup_report.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
