from math import comb
import sympy as s
u,v=s.symbols('u v');M=s.Matrix(3,3,lambda i,j:s.Symbol(f'a{min(i,j)}{max(i,j)}'));E=s.zeros(4);E[1:4,1:4]=M
T2=s.Matrix([[s.Rational(comb(j,i),comb(2,i))if i<=j else 0 for i in range(3)]for j in range(3)])
T3=s.Matrix([[s.Rational(comb(j,i),comb(3,i))if i<=j else 0 for i in range(4)]for j in range(4)])
old=T2*M*T2.T;new=T3*E*T3.T
for i in range(4):
 for j in range(4):
  expected=s.Rational(i*j,9)*old[i-1,j-1]if i*j else 0
  assert s.expand(new[i,j]-expected)==0
assert s.expand((s.Matrix([[1,u,u**2,u**3]])*E*s.Matrix([1,v,v**2,v**3]))[0]-u*v*(s.Matrix([[1,u,u*u]])*M*s.Matrix([1,v,v*v]))[0])==0
print('PASS exact rank-three embedding and ten cubic Bernstein control identities; numerical PSD margin remains an additional restriction')
