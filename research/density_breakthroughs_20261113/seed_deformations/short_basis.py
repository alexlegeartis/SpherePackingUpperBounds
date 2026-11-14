import numpy as np,pathlib,time,sympy as sp
from fpylll import IntegerMatrix,BKZ,LLL
p=pathlib.Path(__file__).resolve().parent/'dim48mod3';D=np.load(p/'reduced.npz');G=D['H'];V=np.load(p/'contacts.npy');n=len(G);rng=np.random.default_rng(1);rng.shuffle(V);W=np.vstack([np.eye(n,dtype=np.int64),V[:200]])
B0=np.rint(np.linalg.cholesky(G)*1e8).astype(np.int64);B=IntegerMatrix.from_matrix((W@B0).tolist());U=IntegerMatrix.identity(len(W));LLL.reduction(B,delta=.99,U=U);U=np.array([[int(U[i,j])for j in range(len(W))]for i in range(len(W))],np.int64);T=(U@W)[-n:];H=T@G@T.T
assert abs(sp.Matrix(T).det())==1
B=IntegerMatrix.from_matrix((T@B0).tolist());U=IntegerMatrix.identity(n);BKZ.reduction(B,BKZ.Param(block_size=30,max_loops=2),U=U,float_type='long double');U=np.array([[int(U[i,j])for j in range(n)]for i in range(n)],np.int64);T=U@T;H=T@G@T.T
print('diagonal',H.diagonal(),flush=True);np.savez(p/'short_reduced.npz',H=H,T=T,G=G);np.savetxt(p/'Qshort.txt',H,fmt='%d');np.savetxt(p/'Rshort.txt',np.linalg.cholesky(H).T,fmt='%.17g')
