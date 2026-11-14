"""Floating discovery: seek a metric descent direction preserving known contacts."""
import sys,time,json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog

prefix=sys.argv[1]
data=np.load(prefix+'.npz'); Q=data['Q']; n=len(Q)
V=np.loadtxt(prefix+'_vectors.txt',dtype=np.int64)
V=np.unique(np.array([v if next((int(z)for z in v if z),0)>0 else -v for v in V]),axis=0)
norms=np.einsum('vi,ij,vj->v',V,Q,V)
assert norms.max()-norms.min()<1e-6
i,j=np.triu_indices(n); factors=np.where(i==j,1.,2.)
A=V[:,i]*V[:,j]*factors
c=np.linalg.inv(Q)[i,j]*factors
c=c/np.linalg.norm(c)
rng=np.random.default_rng(20261114)
ids=rng.choice(len(V),min(5000,len(V)),replace=False)
t=time.time()
for iteration in range(10):
 res=linprog(c,A_ub=-A[ids],b_ub=np.zeros(len(ids)),bounds=(-1,1),method='highs',options={'time_limit':35.})
 print(prefix,'iteration',iteration,'constraints',len(ids),'status',res.message,'objective',res.fun,'seconds',time.time()-t,flush=True)
 if not res.success: break
 bad=np.flatnonzero(A@res.x < -1e-7)
 if not len(bad):
  D=np.zeros((n,n));D[i,j]=res.x;D[j,i]=res.x
  np.savez(prefix+'_descent.npz',D=D,Q=Q,V=V)
  Path(prefix+'_eutaxy.json').write_text(json.dumps({'contacts_pairs':len(V),'objective':res.fun,'constraints':len(ids),'minimum_constraint':float((A@res.x).min()),'iterations':iteration+1},indent=2)+'\n')
  print('VALIDATED ALL CONTACTS objective',res.fun,flush=True);break
 ids=np.unique(np.concatenate([ids,rng.choice(bad,min(3000,len(bad)),replace=False)]))
