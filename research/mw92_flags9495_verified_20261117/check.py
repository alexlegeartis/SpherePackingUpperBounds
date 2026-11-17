from pathlib import Path
from fractions import Fraction as F
from math import factorial,prod,isqrt
import json,hashlib,sys
sys.set_int_max_str_digits(0)
P=Path(__file__).resolve().parent;r=json.load(open(P/'input/results.json'));n=92
# Independently expand Dirichlet(1,...,1) moments of sum(T_i^2).
mu=F(n*factorial(2),n*(n+1))
m2=F(n*factorial(4)+n*(n-1)*factorial(2)**2,prod(range(n,n+4)))
m3=F(n*factorial(6)+3*n*(n-1)*factorial(4)*factorial(2)+n*(n-1)*(n-2)*factorial(2)**3,prod(range(n,n+6)))
var=m2-mu**2;third=m3-3*mu*m2+2*mu**3;a=F(n,2)
C=(n+1)*F(n,n+3)**46*(1+a*(a+1)*var/(2*(1+mu)**2)-a*(a+1)*(a+2)*third/(6*(1+mu)**3))
assert C==F(r['covering_constant']) and C>(n+1)*F(n,n+3)**46
# Alternating-series lower/upper Machin enclosures; odd length is upper.
def atan(q,N):return sum((F((-1)**k,(2*k+1)*q**(2*k+1))for k in range(N)),F(0))
lo=16*atan(5,100)-4*atan(239,101);hi=16*atan(5,101)-4*atan(239,100)
assert 3<lo<hi<F(22,7) and hi==F(r['pi_upper'])
d={int(j):F(v)for j,v in r['depths'].items()}
for j in d:assert 0<d[j]<16 and j*hi**46/factorial(46)*d[j]**46<C*47**7
h1=16-d[1];h2=16-d[3];G=[[F(v)for v in row]for row in r['gram']]
assert G==[[h1,-h1/2],[-h1/2,h2]] and h2>h1>4
D=h1*h2-h1*h1/4;assert D==F(r['determinant'])>0
# G−2I positive definite. Hence any omitted coordinate |x|>=4 or |y|>=4
# has norm>32, larger than every claimed parity minimum.
assert h1>2 and (h1-2)*(h2-2)-h1*h1/4>0
mins={}
for x in range(-3,4):
 for y in range(-3,4):
  if (x,y)==(0,0):continue
  parity=(x%2,y%2);v=h1*x*x-h1*x*y+h2*y*y
  mins[parity]=min(mins.get(parity,v),v)
assert mins=={(0,0):4*h1,(1,0):h1,(0,1):h2,(1,1):h2} and 32>max(mins.values())
# All16 color-pair cases, including wraparound sign changes in horizontal differences.
colors=[(0,0),(1,0),(0,1),(1,1)]
for x in colors:
 for y in colors:
  if x==y:assert mins[(0,0)]>16;continue
  diff=(x[0]-y[0],x[1]-y[1]);par=(diff[0]%2,diff[1]%2)
  depth=d[1] if diff[1]==0 else d[3]
  assert depth+mins[par]>=16
center2=F(2**188,47**14)/D;assert center2==F(r['center_density_squared'])
# Exact outward rounding; no mpmath/float arithmetic enters certification.
safe=F(784434540405174,10**35)
assert safe**2<center2*lo**94/factorial(47)**2
baseline=json.load(open(P/'baseline.json'));old2=F(baseline['row']['density_squared']);assert center2>old2
# Certified gain bracket versus actual preceding public rational certificate.
glo=F(173020,10**5);ghi=F(173021,10**5)
assert (1+glo/100)**2<center2/old2<(1+ghi/100)**2
external2=F(2**184,47**14)/12
ratio2=center2/external2;rr=F(isqrt(ratio2.numerator*10**24//ratio2.denominator),10**12)
assert rr**2<=ratio2<(rr+F(1,10**12))**2
external_gain=[str(100*(rr-1)),str(100*(rr+F(1,10**12)-1))]
report=dict(external_derived_comparator_center_squared=str(external2),external_gain_percent_interval=external_gain,status='PASS exact covering depths, complete vertical parity minima, all color cases, density and rounding; published parent/covering premises retained',dimension=94,safe_occupied_lower=str(safe),safe_decimal='7.84434540405174e-21',determinant=str(D),center_density_squared=str(center2),parity_minima={str(k):str(v)for k,v in mins.items()},gain_over_public_percent_interval=[str(glo),str(ghi)],baseline_commit=baseline['commit'],input_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest()for f in sorted((P/'input').iterdir())})
(P/'RESULT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
