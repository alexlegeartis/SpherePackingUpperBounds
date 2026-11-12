from pathlib import Path
import numpy as np,sympy as sp
p=Path(__file__).resolve().parent;D=np.load(p/'data.npz');G,H,T=[D[k]for k in['G','H','T']]
assert sp.Matrix(G).det()==3**18
assert abs(sp.Matrix(T).det())==1
assert np.array_equal(H,T@G@T.T)
L,B=sp.Matrix(H).LDLdecomposition(hermitian=False);assert all(B[i,i]>0 for i in range(36))
V=np.load(p/'min6.npy');assert V.shape==(576,36);assert np.all(np.einsum('ij,jk,ik->i',V,H,V)==6)
A=V.T@V@H;assert np.array_equal(A@A@A,288*A@A-18432*A)
assert np.trace(A)==3456 and np.trace(A@A)==552960
P=np.load(p/'parts.npy');assert np.array_equal(P.sum(0),18432*H)
assert all(np.array_equal(M,M.T)for M in P)
K=np.rint(np.load(p/'candidate.npy')).astype(np.int64);assert np.array_equal(55296*K,P[0]+3*P[1]+9*P[2])
W=np.load(p/'newshort.npy');assert np.all(np.einsum('ij,jk,ik->i',W,K,W)==6)
assert np.all(np.einsum('ij,jk,ik->i',W,H,W)==14)
print('Exact determinant, basis, positivity, frame polynomial/multiplicities, and rejection witnesses: PASS')
