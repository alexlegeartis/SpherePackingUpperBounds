from pathlib import Path
import numpy as np,sympy as sp,json,time,math
from fpylll import IntegerMatrix,GSO,LLL,BKZ
P=Path(__file__).resolve().parent;source=P.parent/'dual64/source.npz';d=np.load(source);G=d['G'];H=d['H'];Q=IntegerMatrix.from_matrix(G.tolist());U=IntegerMatrix.identity(64);M=GSO.Mat(Q,gram=True,U=U,update=True);lll=LLL.Reduction(M);lll();start=time.time();BKZ.Reduction(M,lll,BKZ.Param(block_size=24,max_loops=3))();T=np.array([[U[i,j] for j in range(64)]for i in range(64)],dtype=np.int64);R=np.array([[Q[i,j]for j in range(64)]for i in range(64)],dtype=np.int64);assert np.array_equal(T@G@T.T,R);Ti=np.array(sp.Matrix(T.tolist()).inv().tolist(),dtype=np.int64);assert np.array_equal(T@Ti,np.eye(64,dtype=np.int64));np.savez_compressed(P/'ne64_primal_reduced.npz',G=G,H=H,T=T,Ti=Ti,R=R);np.savetxt(P/'ne64_primal_Q.txt',R,fmt='%d');r=np.array([M.get_r(i,i)for i in range(64)]);print('BKZ',time.time()-start,'diag',np.diag(R).tolist(),'gso',r.tolist(),flush=True)
for j in range(4,65,4):
 logv=j/2*math.log(math.pi*60)-math.lgamma(j/2+1)-j*math.log(3)-.5*sum(np.log(r[-j:]));print(j,round(logv/math.log(10),3))
