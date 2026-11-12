import numpy as np,sympy as s
from pathlib import Path
p=Path('research/certify_ct36_minimum');f=np.load(p/'bkz36.npz');Q=s.Matrix(f['Q'].tolist());U=s.Matrix(f['U'].tolist());a=np.load('research/ct36_contact_walk_20261112/flex.npz');B=s.Matrix((5*a['G4']+4*a['D']).tolist());assert s.Matrix(np.load('research/ct36_contact_walk_20261112/candidate_1.npy').tolist())==B
assert all(x.is_Integer for x in Q) and Q==Q.T and abs(U.det())==1 and Q==U*B*U.T
assert np.array_equal(np.loadtxt(p/'Q36.txt',dtype=np.int64),np.array(Q).astype(np.int64))
L,D=Q.LDLdecomposition(hermitian=False);rows=np.loadtxt(p/'ldl_bounds.txt',dtype=np.int64);ix=0
for i in range(36):
 for x in [D[i,i]]+[L[j,i] for j in range(36)]:
  lo,hi=map(int,rows[ix]);ix+=1;assert s.Rational(lo,2**40)<=x<=s.Rational(hi,2**40);assert max(abs(lo),abs(hi))<2**53
 assert D[i,i]>0
print('PASS exact candidate=20G(1/5), Q36 text identity, U unimodularity, rational interval enclosures')
inv=Q.inv();print('Max genuine radius159 vector coordinate bound',max(float((159*inv[i,i])**s.Rational(1,2)) for i in range(36)));print('Max Q',max(abs(x) for x in Q))
