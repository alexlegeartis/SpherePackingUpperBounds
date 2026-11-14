"""Numerical eutaxy screen using explicit, exact-length section contacts.

A positive result is only a discovery observation until rationally certified.
"""
from pathlib import Path
import numpy as np,json,time
from scipy.optimize import nnls
from scipy.linalg import lstsq
P=Path(__file__).resolve().parent;ROOT=P.parents[2]
F=np.load(P/'full_contacts.npz');V=F['V'];types=F['oldnorm']
Q=np.loadtxt(ROOT/'constructions/dim36/data/Q_original.txt',dtype=np.int64)
out=[]
for idx in [20,29,30,33,34,35]:
 keep=[i for i in range(36)if i!=idx];mask=V[:,idx]==0;W=V[mask][:,keep];typ=types[mask]
 G=Q[np.ix_(keep,keep)]/20;C=np.linalg.cholesky(G);X=W@C/np.sqrt(8)
 i,j=np.triu_indices(35);mul=np.where(i==j,1,np.sqrt(2));A=X[:,i]*X[:,j]*mul
 y=(np.eye(35)[i,j]*mul)
 levels=np.unique(typ);T=np.stack([A[typ==v].sum(axis=0)for v in levels],axis=1)
 weights0,res0=nnls(T,y);w=np.zeros(len(W))
 for v,a in zip(levels,weights0):w[typ==v]=a
 H=A.T@A;r=y-A.T@w;z,res,rank,sv=lstsq(H,r,cond=1e-11);w+=A@z
 err=float(np.max(np.abs(A.T@w-y)))
 row={'deleted_coordinate':idx,'contacts_up_to_sign':len(W),'constraint_rank':int(rank),'minimum_weight':float(w.min()),'maximum_weight':float(w.max()),'equality_error':err,'condition':float(sv[0]/sv[-1])}
 print(json.dumps(row),flush=True);out.append(row)
 np.savez(P/f'eutaxy_{idx}.npz',weights=w,W=W,Q=Q[np.ix_(keep,keep)])
 (P/'eutaxy_screen.json').write_text(json.dumps(out,indent=2))
