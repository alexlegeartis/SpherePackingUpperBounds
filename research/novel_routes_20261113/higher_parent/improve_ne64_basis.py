from pathlib import Path
import numpy as np,sympy as sp,math,time,json
from fpylll import IntegerMatrix,GSO,LLL,BKZ
P=Path(__file__).resolve().parent;d=np.load(P/'ne64_primal_reduced.npz');G=d['G'];H=d['H'];rng=np.random.default_rng(266333);best=1e99;hist=[]
for seed in range(15):
 S=np.eye(64,dtype=np.int64)[rng.permutation(64)];Q=IntegerMatrix.from_matrix((S@G@S.T).tolist());U=IntegerMatrix.identity(64);M=GSO.Mat(Q,gram=True,U=U,update=True);lll=LLL.Reduction(M);lll();t=time.time();BKZ.Reduction(M,lll,BKZ.Param(block_size=28 if seed<12 else 32,max_loops=3))();T=np.array([[U[i,j]for j in range(64)]for i in range(64)],dtype=np.int64)@S;R=T@G@T.T;r=np.array([M.get_r(i,i)for i in range(64)]);logs=[j/2*math.log(math.pi*60)-math.lgamma(j/2+1)-j*math.log(3)-.5*sum(np.log(r[-j:])) for j in range(1,65)];cost=math.log(sum(math.exp(x-max(logs))for x in logs))+max(logs);row={'seed':seed,'estimated_log10_total':cost/math.log(10),'seconds':time.time()-t};hist.append(row);print(row,flush=True)
 if cost<best:
  best=cost;Ti=np.array(sp.Matrix(T.tolist()).inv().tolist(),dtype=np.int64);assert np.array_equal(T@Ti,np.eye(64,dtype=np.int64));np.savez_compressed(P/'ne64_primal_best.npz',G=G,H=H,T=T,Ti=Ti,R=R);np.savetxt(P/'ne64_primal_best_Q.txt',R,fmt='%d')
(P/'basis_search.json').write_text(json.dumps(hist,indent=2))
