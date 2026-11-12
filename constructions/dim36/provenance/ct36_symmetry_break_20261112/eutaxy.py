import numpy as np, scipy.optimize as so, json, time
from pathlib import Path
P=Path(__file__).resolve().parent
N=np.load(P/'reduced.npz')['N'];H=N/4;V=np.loadtxt(P/'vectors.txt',dtype=np.int64).reshape(-1,36);V=np.unique(np.array([v if next(x for x in v if x)!=abs(next(x for x in v if x)) else -v for v in V]),axis=0)
rng=np.random.default_rng(123);V=V[rng.choice(len(V),min(4000,len(V)),replace=False)]
i,j=np.triu_indices(36);A=(V[:,i]*V[:,j]).T.astype(float);b=np.linalg.inv(H)[i,j]
t=time.time();r=so.linprog(np.zeros(len(V)),A_eq=A,b_eq=b,bounds=(0,None),method='highs',options={'time_limit':25})
print(r.message,'elapsed',time.time()-t,flush=True)
if r.success:
 ids=np.where(r.x>1e-9)[0];print('support',len(ids),'residual',np.max(abs(A@r.x-b)))
 np.savez(P/'eutaxy_candidate.npz',vectors=V[ids],weights=r.x[ids],N=N)
json.dump({'message':r.message,'success':bool(r.success),'tested_rows':len(V),'method':'floating cone feasibility, not exact proof'},open(P/'eutaxy_result.json','w'),indent=2)
