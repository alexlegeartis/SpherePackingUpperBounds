import json,sympy as s
from fractions import Fraction as F
from pathlib import Path
p=Path(__file__).parent;x=s.symbols('x');polys={}
for k,f in [(37,'rank2_candidate37_B17677213_20261117.json')]:
 co=json.loads((p/f).read_text())['even_coefficients'];a=s.Rational(k-2,2);polys[k]=s.Poly(sum(s.Rational(str(v))*s.gegenbauer(2*j,a,x)/s.gegenbauer(2*j,a,1)for j,v in enumerate(co))-s.Rational(2,10000),x)
rows=[]
for n,k in [(37,37)]:
 Q=polys[k];rest=Q;co=[]
 for j in range(Q.degree()//2,-1,-1):
  basis=s.Poly(s.gegenbauer(2*j,s.Rational(n-2,2),x)/s.gegenbauer(2*j,s.Rational(n-2,2),1),x);a=rest.nth(2*j)/basis.nth(2*j);assert a>=0;co.append(a);rest-=a*basis
 assert rest.is_zero and co[-1]>0;qn=F(Q.eval(1)/co[-1]);old=json.loads((p/f'certificate{n}.json').read_text());a=list(map(F,old['cyclic_angular']['gegenbauer_coefficients']));B=sum(a);assert B==F(old['cyclic_angular']['B_exact']);qold=sum(a[::2]);C=F(old['C']);r=F(old.get('origin_ratio','1'));origin0=F(old.get('origin_C',old['C']));assert origin0==C*r*r
 # The published38/39 radial verifier explicitly uses origin C and cross -2C/B, hence r=1.
 if n in (38,39):assert 'origin_ratio' not in old and r==1
 c=C/B**2;u0=c*qold;v0=c*r*B;Gamma=(1-F(1,100))*(1-F(1,10000));u=c*B*Gamma/2 if n in (38,39) else u0;v=v0+u-u0
 assert u>=u0>0 and 2*u<=c*B*Gamma
 origin=qn*v*v/u;extra=2*(2**n-1)*F(old['kappa'])+F(old['extra_objective']);obj=F(old['p'][0])+origin0+extra;newobj=F(old['p'][0])+origin+extra;upper=F(old['safe_density'])*newobj/obj;m=10**(14 if n==38 else 15);safe=F((upper.numerator*m+upper.denominator-1)//upper.denominator,m)
 assert upper<=safe<F(old['safe_density']);row=dict(n=n,kernel=k,safe=str(safe),safe_decimal=format(float(safe),'.15f'),decrease_percent=float(100*(1-newobj/obj)));rows.append(row);print(row)
(p/'PROSPECTIVE_RESULTS_37.json').write_text(json.dumps(rows,indent=2));print('Exact transfer, PSD and rounding PASS; independent gate/dependency audit still required.')
