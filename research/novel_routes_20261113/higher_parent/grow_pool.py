import numpy as np,json,time
from pathlib import Path
p=Path(__file__).resolve().parent;G=np.load(p/'G72.npy');E=np.eye(72,dtype=np.int16)
V=[x for x in E]
for i in range(72):
 for j in range(i+1,72):
  if abs(G[i,j])==4:V.append(E[i]-int(np.sign(G[i,j]))*E[j])
V=np.array(V);V*=np.sign(V[np.arange(len(V)),np.argmax(V!=0,axis=1)])[:,None];V=np.unique(V,axis=0);rng=np.random.default_rng(713267)
for step in range(10):
 N=len(V);sel=rng.choice(N,min(320,N),replace=False);C=(V[sel].astype(np.int64)@G@V.T)
 i,j=np.where(np.abs(C)==4);W=V[sel[i]]-np.sign(C[i,j])[:,None]*V[j];W=W.astype(np.int16);W*=np.sign(W[np.arange(len(W)),np.argmax(W!=0,axis=1)])[:,None];V=np.unique(np.vstack([V,W]),axis=0)
 if len(V)>50000:V=V[rng.choice(len(V),50000,replace=False)]
 print(step,len(V),flush=True)
 if len(V)>=50000:break
assert np.all(np.einsum('bi,ij,bj->b',V.astype(np.int64),G,V.astype(np.int64))==8)
np.savez(p/'gamma72_norm8_pool.npz',G=G,V=V)
