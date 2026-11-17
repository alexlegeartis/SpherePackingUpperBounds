from pathlib import Path
from fractions import Fraction as F
from math import factorial
from functools import lru_cache
from flint import fmpq as Q,fmpq_poly as Poly
import json,sys,hashlib
if not __debug__:raise RuntimeError('Assertions required')
sys.set_int_max_str_digits(0);D=Path(__file__).resolve().parent
atan=lambda q,N:sum((F((-1)**j,(2*j+1)*q**(2*j+1))for j in range(N)),F(0));piL=16*atan(5,200)-4*atan(239,201)
reports=[]
for n in [49]:
 P=D/f'dim{n}';r=json.load(open(P/'repaired.json'));c=json.load(open(P/'annular_final_certificate.json'));ang=json.load(open(P/'angular.json'));p=Poly(list(map(Q,r['P'])));q=Poly(list(map(Q,r['Q'])))
 assert r['dimension']==n and F(r['contact_x'])==22 and p.degree()==q.degree()==71
 # Independent Laguerre recurrence rather than producer's monomial/rising-product transform.
 x=Poly([0,1]);alpha=Q(n,2)-1;L=[Poly([1]),Poly([alpha+1,-1])]
 for k in range(2,72):L.append(((2*k-1+alpha-x)*L[-1]-(k-1+alpha)*L[-2])/k)
 assert sum((p[k]*factorial(k)*L[k]for k in range(72)),Poly([0]))==q and q(0)==1
 fac=Poly([1])
 for f in r['factors']:
  if f['type']=='negative_linear':
   a=Q(f['a']);assert a>0;v=Poly([1,1/a])
  elif f['type']=='positive_double':
   a=Q(f['mid']);assert a>0;v=Poly([1,-1/a])**2
  elif f['type']=='complex_pair':
   re,im=Q(f['re']),Q(f['im']);assert im>0;v=(Poly([-re,1])**2+im*im)/(re*re+im*im)
  else:raise AssertionError(f['type'])
  fac*=v
 assert fac==q
 co=list(map(Q,ang['gegenbauer_coefficients']));assert co[0]==1 and min(co)>=0;B=sum(co);g=[Poly([1]),x]
 for k in range(2,len(co)):g.append(((2*k+n-4)*x*g[-1]-(k-1)*g[-2])/(k+n-3))
 assert sum((a*b for a,b in zip(co,g)),Poly([0]))==Poly(list(map(Q,ang['polynomial_coefficients'])))
 assert B==Q(c['q']) and F(ang['profile_a_exact'])==-192 and F(ang['profile_b_exact'])==216
 budget=F(0);last=F(22);assert len(c['rows'])==632
 @lru_cache(None)
 def ratio(k):return F(factorial(n-2+k)*factorial((n-3)//2),factorial(n-2)*factorial((n-3)//2+k))
 for row in c['rows']:
  l,r=F(row['l']),F(row['r']);assert l==last and r-l==F(1,4);last=r;s=1-F(11)/r;k=row['cap_degree'];assert type(k)==int and k>=0
  gap=F(k,k+n-1)-s;assert gap>0;bound=(1-s)*ratio(k)/gap;assert row['cap']==bound.numerator//bound.denominator
  h=F(row['height']);assert h>=0;budget+=h*row['cap']
 assert last==180 and budget==F(c['budget'])
 a,z,s=[F(c[k])for k in ['a','z','s']];assert a>0 and s==z*z/a+F(1,10**12) and a*s-z*z>0 and a+2*z<0
 assert c['P']==list(map(str,p));objective=F(str(p[0]))+s+budget;assert objective>0
 num=4*objective**2*22**n*F(factorial((n+1)//2),factorial(n+1))**2
 safe=F({49:'2.87119915e-7',51:'1.54978656e-7'}[n]);assert num<safe**2*piL
 reports.append(dict(dimension=n,status='PASS',strict_upper=str(safe),safe_decimal=float(safe),annuli=632,nodes=c['nodes'],budget=float(budget),angular_sha256=hashlib.sha256((P/'angular.json').read_bytes()).hexdigest(),odd_normalization='independent exact square and200/201-term Machin lower pi'))
(D/'exact_result.json').write_text(json.dumps(reports,indent=2));print(json.dumps(reports,indent=2))
