"""Independent centered Taylor enclosure using repeated differential operators."""
from fractions import Fraction as F
from math import comb,factorial
from flint import fmpq as Q,arb,ctx
ctx.prec=512

def ar(q):return arb(q.numerator)/q.denominator
def shift(p,a):
 out=[p[-1]]
 for v in reversed(p[:-1]):out=[a*out[0]+v]+[out[j-1]+a*out[j] for j in range(1,len(out))]+[out[-1]]
 return out

def upper(p,tau,C,B,qminus,power,b,l,u,origin_ratio):
 mid=(l+u)/2;r=(u-l)/2;Z=lambda v:Q(v.numerator,v.denominator)
 m=Z(mid);rr=Z(r);tau=Z(tau);C=Z(C);B=Z(B);qm=Z(qminus);b=Z(b)
 rad=shift([Z(v)*tau**j for j,v in enumerate(p)],m)
 specs=[(Q(1),tau/2,Q(0),rad)]
 for scale,beta,const,degree in [(C/B*(1+qm/B),2*b,2*b,2*power),(-2*Z(origin_ratio)*C/B,b,b,power),(2*C/B*4**power,5*b,2*b,2*power)]:
  specs.append((scale,beta,const,[Q(comb(degree,j))*m**(degree-j) for j in range(degree+1)]))
 co=[arb(0) for _ in range(17)];rem=arb(0)
 for scale,beta,const,poly in specs:
  amp=ar(scale)*(ar(const-beta*m)).exp()
  for j in range(18):
   if j<17:co[j]+=amp*ar(poly[0]/factorial(j))
   else:
    bound=sum((abs(v)*rr**k for k,v in enumerate(poly)),Q(0))
    rem+=ar(abs(scale))*(ar(const-beta*Z(l))).exp()*ar(bound*rr**17/factorial(17))
    break
   poly=[(k+1)*poly[k+1]-beta*poly[k] for k in range(len(poly)-1)]+[-beta*poly[-1]]
 # Taylor in s=t-mid; shift to t-l and convert to Bernstein on[0,2r].
 shifted=[sum((co[k]*ar(Q(comb(k,j))*(-rr)**(k-j)) for k in range(j,17)),arb(0)) for j in range(17)]
 return [sum((shifted[j]*ar((2*rr)**j*Q(comb(i,j),comb(16,j))) for j in range(i+1)),arb(0))+rem for i in range(17)]
