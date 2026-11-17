"""Independent exact LP-polynomial, volume and section verification.
No Jacobi generator/norm formula is needed: verify the supplied certificate itself.
"""
from pathlib import Path
from fractions import Fraction as F
from math import factorial,prod
import json,sympy as S
P=Path(__file__).resolve().parent;x=S.symbols('x')
r=json.loads((P/'kl_report.json').read_text());gammas={}
# Independently use 120/121 terms rather than producer100/101.
def arctan(q,N):return sum((F((-1)**k,(2*k+1)*q**(2*k+1)) for k in range(N)),F(0))
piL=16*arctan(5,120)-4*arctan(239,121)
for a in r['hermite_bounds']:
 n=a['dimension'];s=S.Rational(a['s']);co=list(map(S.Rational,a['nonnegative_normalized_gegenbauer_coefficients']))
 assert n in (90,91,92) and len(co)==12 and co[0]>0 and all(c>=0 for c in co) and -1<s<=S.Rational(1,2)
 # Normalized Gegenbauer recurrence built separately from special-function calls.
 Q=[S.Poly(1,x,domain=S.QQ),S.Poly(x,x,domain=S.QQ)]
 for j in range(1,11):Q.append(((2*j+n-2)*S.Poly(x,x)*Q[j]-j*Q[j-1]).mul_ground(S.Rational(1,j+n-2)))
 f=sum((c*q for c,q in zip(co,Q)),S.Poly(0,x))
 quotient,remainder=S.div(f,S.Poly(x-s,x));assert remainder.is_zero
 # f=(x-s)*positive rational*polynomial square, hence nonpositive for x<=s.
 factor,terms=S.factor_list(quotient);assert factor>0 and all(e%2==0 for _,e in terms)
 A=F(str(f.eval(1)/co[0]));assert A==F(a['A_bound'])
 U2=A*A*((1-F(str(s)))/2)**n;assert U2==F(a['sphere_upper_squared'])
 V2=piL**n/F(factorial(n//2)**2) if n%2==0 else F(2**(n+1),prod(range(1,n+1,2))**2)*piL**(n-1)
 gamma=F(a['hermite_upper']);assert gamma**n>=4**n*U2/V2;gammas[n]=gamma
D=F(2**184,47**14)
for a in r['sections']:
 n=a['dimension']+1;floor=F(a['center_density_squared_lower']);assert floor.denominator==1
 assert floor**n*gammas[n]**n<=D**(n-1)<(floor+1)**n*gammas[n]**n
 D=floor
print('PASS independent recurrence, polynomial square/sign, PSD coefficients, pi enclosure, Hermite and section powers')
