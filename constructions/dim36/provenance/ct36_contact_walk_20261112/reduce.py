import numpy as np,sympy as sp,time
from pathlib import Path
P=Path(__file__).resolve().parent
for sg in[1,-1]:
 Q=np.load(P/f'candidate_{sg}.npy');B=np.rint(np.linalg.cholesky(Q/20)*10**7).astype(np.int64);t=time.time();_,U=sp.Matrix(B).lll_transform(delta=sp.Rational(99,100));U=np.array(U,dtype=np.int64);R=U@Q@U.T;np.savez(P/f'reduced_{sg}.npz',Q=R,U=U);np.savetxt(P/f'Qred{sg}.txt',R,fmt='%d');np.savetxt(P/f'Rred{sg}.txt',np.linalg.cholesky(R/20).T);print('LLL',sg,'seconds',time.time()-t,'diag',R.diagonal().min(),R.diagonal().max(),flush=True)
