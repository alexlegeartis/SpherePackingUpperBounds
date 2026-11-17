from pathlib import Path
from fractions import Fraction as F
from math import factorial,comb
from functools import lru_cache
import json,sys,hashlib
from flint import fmpq as Q,fmpq_poly as Poly
if not __debug__:raise RuntimeError('Assertions required')
sys.set_int_max_str_digits(0);D=Path(__file__).resolve().parent
r=json.load(open(D/'repaired.json'));c=json.load(open(D/'annular_final_certificate.json'));ang=json.load(open(D/'angular.json'));p=Poly(list(map(Q,r['P'])));q=Poly(list(map(Q,r['Q'])))
assert r['dimension']==50 and F(r['contact_x'])==22
# Independently use the direct binomial monomial transform.
assert q==Poly([sum((p[k]*(-1)**j*factorial(k)*comb(k+24,k-j)/factorial(j)for k in range(j,72)),Q(0))for j in range(72)]) and q(0)==1
factor=Poly([1])
for row in r['factors']:
 t=row['type']
 if t=='negative_linear':
  x=Q(row['a']);assert x>0;f=Poly([1,1/x])
 elif t=='positive_double':
  x=Q(row['mid']);assert x>0;f=Poly([1,-1/x])**2
 elif t=='complex_pair':
  re,im=Q(row['re']),Q(row['im']);assert im>0;f=(Poly([-re,1])**2+im*im)/(re*re+im*im)
 else:raise AssertionError(t)
 factor*=f
assert factor==q
# Exact normalized Gegenbauer reconstruction via three-term recurrence.
x=Poly([0,1]);g=[Poly([1]),x]
co=list(map(Q,ang['gegenbauer_coefficients']));assert co[0]==1 and min(co)>=0
for k in range(2,len(co)):g.append(((2*k+46)*x*g[-1]-(k-1)*g[-2])/(k+47))
qa=sum((a*b for a,b in zip(co,g)),Poly([0]));assert qa==Poly(list(map(Q,ang['polynomial_coefficients'])))
assert sum(co)==qa(1)==Q(c['q']) and F(ang['profile_a_exact'])==-192 and F(ang['profile_b_exact'])==216
# Closed factorial formula, independent of producer's rising-product cap calculation.
@lru_cache(None)
def moment_ratio(k):return F(4**k*factorial(48+k)*factorial(24+k),factorial(48+2*k)*factorial(24))
assert len(c['rows'])==632;budget=F(0);last=F(22);maxk=0
for row in c['rows']:
 l,r=F(row['l']),F(row['r']);assert l==last and r-l==F(1,4);last=r
 s=1-F(11)/r;k=row['cap_degree'];assert type(k)==int and k>=0
 gap=F(k,k+49)-s;assert gap>0
 bound=(1-s)*moment_ratio(k)/gap;assert row['cap']==bound.numerator//bound.denominator
 h=F(row['height']);assert h>=0;budget+=row['cap']*h;maxk=max(maxk,k)
assert last==180 and budget==F(c['budget'])
a,z,s=[F(c[k])for k in ['a','z','s']];assert a>0 and s==z*z/a+F(1,10**12) and a*s>z*z and a+2*z<0
assert c['P']==list(map(str,p))
density=F(22,4)**25/F(factorial(25))*(F(str(p[0]))+s+budget)
assert 0<density<F('2.06573053e-7')
report=dict(status='PASS',dimension=50,scope='unrestricted sphere packings',strict_upper='2.06573053e-7',density_decimal=float(density),budget_decimal=float(budget),annuli=632,nodes=c['nodes'],max_cap_degree=maxk,angular_sha256=hashlib.sha256((D/'angular.json').read_bytes()).hexdigest())
(D/'independent_exact_result.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
