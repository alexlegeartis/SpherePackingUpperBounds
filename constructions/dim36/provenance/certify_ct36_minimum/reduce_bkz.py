import sys,time
sys.path.insert(0,'/tmp/ct36_svp_deps')
from fpylll import IntegerMatrix, BKZ, LLL
import numpy as np,sympy as sp
from pathlib import Path
p=Path(__file__).parent;q=np.load(p.parent/'ct36_contact_walk_20261112/candidate_1.npy');b=IntegerMatrix.from_matrix(np.rint(np.linalg.cholesky(q)*1e8).astype(np.int64).tolist());u=IntegerMatrix.identity(36)
for block in [20,28,36]:
 t=time.time();BKZ.reduction(b,BKZ.Param(block_size=block,max_loops=3),U=u)
 U=np.array([[int(u[i,j])for j in range(36)]for i in range(36)],dtype=np.int64);Q=U@q@U.T
 print(block,time.time()-t,'min',Q.diagonal().min(),'diag',Q.diagonal(),flush=True)
 np.savez(p/f'bkz{block}.npz',U=U,Q=Q);np.savetxt(p/f'Q{block}.txt',Q,fmt='%d');np.savetxt(p/f'R{block}.txt',np.linalg.cholesky(Q/20).T)
 # reset B/U for next reduction so transforms do not lose provenance
 b=IntegerMatrix.from_matrix(np.rint(np.linalg.cholesky(q)*1e8).astype(np.int64).tolist());u=IntegerMatrix.identity(36)
