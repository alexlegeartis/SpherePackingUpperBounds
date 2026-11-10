"""Independent origin-coordinate tightening audit. Inherits source sign proofs."""
from fractions import Fraction as F
from pathlib import Path
import json
from flint import arb, ctx
assert __debug__
ctx.prec=256
P=Path(__file__).resolve().parent
raw=json.loads((P/'witness.json').read_text(),parse_float=str)
Q=[[F(x) for x in row] for row in json.loads((P/'witness_psd_exact.json').read_text())['matrices'][0]]
assert raw['dimension']==5 and raw['J']==1 and raw['p']==4 and len(Q)==6
# Fraction-only LDL and triangular solves, independently of producer's SymPy inverse.
B=[r[1:] for r in Q[1:]];b=[r[0] for r in Q[1:]];n=len(B)
L=[[F(i==j) for j in range(n)] for i in range(n)];d=[]
for j in range(n):
 d.append(B[j][j]-sum(L[j][k]**2*d[k] for k in range(j)))
 assert d[-1]>0
 for i in range(j+1,n):
  L[i][j]=(B[i][j]-sum(L[i][k]*L[j][k]*d[k] for k in range(j)))/d[j]
assert all(B[i][j]==sum(L[i][k]*d[k]*L[j][k] for k in range(n)) for i in range(n) for j in range(n))
y=[]
for i in range(n): y.append(b[i]-sum(L[i][j]*y[j] for j in range(i)))
z=[y[i]/d[i] for i in range(n)];x=[F(0)]*n
for i in reversed(range(n)):x[i]=z[i]-sum(L[j][i]*x[j] for j in range(i+1,n))
assert all(sum(B[i][j]*x[j] for j in range(n))==b[i] for i in range(n))
alpha=Q[0][0]-sum(b[i]*x[i] for i in range(n));assert alpha>0
W=[x]+[[F(i==j) for j in range(n)]for i in range(n)]
assert all(Q[i][j]-(alpha if i==j==0 else 0)==sum(W[i][a]*B[a][b]*W[j][b] for a in range(n) for b in range(n)) for i in range(n+1) for j in range(n+1))
c=list(map(F,raw['x'][:32]));H=sum((-1)**i*v for i,v in enumerate(c));assert H>0
obj=(sum(c)+Q[0][0]+sum(F(a)*F(b) for a,b in zip(raw['count_weights'],raw['count_bounds']))-alpha)/H
def ar(q):return arb(q.numerator)/arb(q.denominator)
density=arb.pi()**2*ar(obj)/60
assert density<arb('0.510501095')
print('PASS: independent Fraction LDL, Schur solves, residual factorization and Arb strict rounding')
print('alpha =',alpha)
print('density =',density)
