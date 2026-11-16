"""Exact symbolic checks for the dimension-uniform determinant proof."""
import sympy as s
n,t,k=s.symbols('n t k');b=4/((n-1)*(n-2))
N=(n-2)-2*(n*n-4*n+5)*t+2*n*(n-3)*t*t
numer=s.expand(-(1-(n-3)*t)*(1+2*t)-(n-3)*(1-t)*(1+2*t)+2*(n-2)*(1-t)*(1-(n-3)*t))
assert s.expand(numer-N)==0
P=n**5-16*n**4+81*n**3-158*n**2+116*n-88
Q=n**4-7*n**3+11*n**2+n+10
assert s.simplify(N.subs(t,b)-P/((n-2)**2*(n-1)**2))==0
assert s.simplify(s.diff(N,t).subs(t,b)+2*Q/((n-2)*(n-1)))==0
for p in [P,Q,n*n-7*n+14]:
 coeff=s.Poly(s.expand(p.subs(n,k+9)),k).all_coeffs()
 assert all(x>0 for x in coeff)
 print('Positive n>=9 via shifted coefficients:',coeff)
for nn in range(9,97):
 bb=s.Rational(4,(nn-1)*(nn-2));assert (nn-3)*bb<1
 assert N.subs({n:nn,t:bb})>0
 ss=s.Rational(8,nn-1)
 assert 1+ss/2-(nn-1)*ss*ss/16==1
print('All symbolic identities and dimension-uniform positivity certificates PASS')
