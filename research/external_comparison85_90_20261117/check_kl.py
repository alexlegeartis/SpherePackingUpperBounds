from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext,ROUND_CEILING
from math import factorial,prod
import sympy as S,json
p=Path(__file__).parent;x=S.symbols('x');out=[]
def atan(q,N):return sum((F((-1)**k,(2*k+1)*q**(2*k+1))for k in range(N)),F(0))
pl=16*atan(5,100)-4*atan(239,101);ph=16*atan(5,101)-4*atan(239,100)
def v2(n,pi):return pi**n/F(factorial(n//2)**2) if n%2==0 else F(2**(n+1),prod(range(1,n+1,2))**2)*pi**(n-1)
for row in json.load(open(p/'search.json')):
 n=row['n'];k=row['k'];assert row['kind']=='odd';s=S.Rational(round(row['s']*10**9),10**9);alpha=S.Rational(n-1,2);beta=S.Rational(n-3,2);K=0
 for i in range(k):
  P=S.jacobi(i,alpha,beta,x);norm=(alpha+beta+1)/(2*i+alpha+beta+1)*S.rf(alpha+1,i)*S.rf(beta+1,i)/(S.factorial(i)*S.rf(alpha+beta+1,i));K+=P*P.subs(x,s)/norm
 f=S.Poly(S.expand((x-s)*K**2),x);rem=f;coeff=[S.Rational(0)]*(f.degree()+1)
 for j in range(f.degree(),-1,-1):
  Q=S.Poly(S.gegenbauer(j,S.Rational(n-2,2),x)/S.gegenbauer(j,S.Rational(n-2,2),1),x);c=rem.nth(j)/Q.nth(j);coeff[j]=c;rem-=c*Q
 assert rem.is_zero and coeff[0]>0 and all(c>=0 for c in coeff);assert -1<s<S.Rational(1,2)
 A=f.eval(1)/coeff[0];U2=F(str(A*A*((1-s)/2)**n));gp=4**n*U2/v2(n,pl)
 with localcontext()as c:
  c.prec=100;g=(Decimal(gp.numerator)/gp.denominator).ln()/n;g=g.exp();ceil=int((g*10**12).to_integral_value(rounding=ROUND_CEILING));gamma=F(ceil,10**12)
 assert gamma**n>=gp
 out.append({'dimension':n,'k':k,'s':str(s),'A_bound':str(A),'sphere_upper_squared':str(U2),'hermite_upper':str(gamma),'nonnegative_normalized_gegenbauer_coefficients':list(map(str,coeff))})
 print('CERTIFIED',n,'gamma <=',float(gamma),flush=True)
D=F(2**184,47**14);sections=[];by={a['dimension']:a for a in out}
for n in [92,91,90]:
 gamma=F(by[n]['hermite_upper']);target=D**(n-1)/gamma**n
 with localcontext()as c:
  c.prec=100;v=((Decimal(target.numerator)/target.denominator).ln()/n).exp();q=int(v)
 while F(q)**n>target:q-=1
 while F(q+1)**n<=target:q+=1
 assert F(q)**n<=target<F(q+1)**n;D=F(q);sections.append({'dimension':n-1,'center_density_squared_lower':str(D),'source_hermite_upper':str(gamma)})
 print('SECTION',n-1,'center square floor',q,flush=True)
(p/'kl_report.json').write_text(json.dumps({'source':'Afkhami-Jeddi et al.2020 Eqs5.10,5.16; Levenshtein spherical LP and Cohn-Zhao transfer','method':'Fixed rational-angle published construction, exact polynomial positivity and arithmetic verification; no claim of global optimization','hermite_bounds':out,'sections':sections},indent=2)+'\n')
