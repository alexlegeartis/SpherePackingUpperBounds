import sympy as s,numpy as np,json
from pathlib import Path
p=Path('research/ct36_contact_walk_20261112/flex.npz'); x=np.load(p);N=s.Matrix(x['G4'].tolist());D=s.Matrix(x['D'].tolist());G=N/4
assert G==G.T and D==D.T
print('detG',G.det(), 'expected',s.Integer(3)**18*s.Rational(5,4)**12,flush=True)
A=N.inv()*D*4
print('trace',s.trace(A),'traceA2',s.trace(A*A),flush=True)
print('A2 scalar identity?',A*A==s.eye(36)*s.Rational(16,15),flush=True)
print('charpoly',s.factor(A.charpoly().as_expr()),flush=True)
V=x['V'];di=np.array(x['D'],dtype=np.int64);ni=np.array(x['G4'],dtype=np.int64)
q=np.einsum('bi,ij,bj->b',V,ni,V);r=np.einsum('bi,ij,bj->b',V,di,V)
print('V shape / q range / r range',V.shape,q.min(),q.max(),r.min(),r.max(),flush=True)
print('rank D',D.rank(),flush=True)
L,B=G.LDLdecomposition(hermitian=False)
assert all(B[i,i]>0 for i in range(36))
assert A**3==s.Rational(16,15)*A
ratio=s.Rational(12288,15625)/(1-s.Rational(16,15)*s.Rational(1,5)**2)**6
print('SPD exact LDL confirmed; cubic confirmed')
print('ratio if min8 at t1/5',ratio,float(ratio),flush=True)
