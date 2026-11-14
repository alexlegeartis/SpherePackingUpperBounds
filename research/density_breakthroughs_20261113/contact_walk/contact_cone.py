"""Cutting-plane search for determinant descent; numerical exploration only."""
import json,time
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
P=Path(__file__).resolve().parent
raw=np.loadtxt(P/'contacts_reduced.txt',dtype=np.int64)
Q=np.loadtxt(P/'Q.txt',dtype=np.int64)
V=raw[:,1:]
assert np.all(raw[:,0]==160)
assert np.all(np.einsum('bi,ij,bj->b',V,Q,V)==160)
n=len(Q);i,j=np.triu_indices(n);mult=np.where(i==j,1,2)
gi=np.linalg.inv(Q);g=gi[i,j]*mult
rng=np.random.default_rng(20261113)
ids=set(map(int,rng.choice(len(V),2500,replace=False)))
history=[]
for it in range(16):
 rows=V[sorted(ids)];A=(rows[:,i]*rows[:,j]*mult).astype(float)
 t=time.time();r=linprog(g,A_ub=-A,b_ub=np.zeros(len(A)),bounds=(-1,1),method='highs',options={'time_limit':50})
 item={'iteration':it,'rows':len(A),'elapsed':time.time()-t,'message':r.message,'objective':None if r.x is None else float(g@r.x)}
 if r.x is None:history.append(item);break
 D=np.zeros_like(Q,dtype=float);D[i,j]=r.x;D[j,i]=r.x
 viol=np.einsum('bi,ij,bj->b',V,D,V)
 item['worst_contact_slope']=float(viol.min());item['negative_contact_count']=int((viol<-1e-7).sum())
 print(json.dumps(item),flush=True);history.append(item)
 np.savez(P/f'cone_{it}.npz',Q=Q,D=D,ids=np.array(sorted(ids)),dual=r.ineqlin.marginals)
 if g@r.x>-1e-8 or viol.min()>-1e-7:break
 bad=np.where(viol<-1e-7)[0];order=bad[np.argsort(viol[bad])]
 ids.update(map(int,order[:500]))
 if len(bad)>500:ids.update(map(int,rng.choice(bad,500,replace=False)))
(P/'cone_history.json').write_text(json.dumps(history,indent=2)+'\n')
