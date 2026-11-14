"""Screen inherited sections for genuine infinitesimal metric flexibility.

Discovery only: contact vectors are explicit but the shell list is incomplete.
Any candidate needs a new complete minimum proof before use as a bound.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import linprog

P=Path(__file__).resolve().parent
ROOT=P.parents[2]
BASE=ROOT/'constructions/dim36'
Q=np.loadtxt(BASE/'data/Q_original.txt',dtype=np.int64)
G=Q/20
f=np.load(BASE/'provenance/ct36_contact_walk_20261112/flex.npz')
V=f['V'].astype(np.int64)
assert np.array_equal(f['G4']+f['D'],Q//5)
assert np.all(np.einsum('bi,ij,bj->b',V,Q,V)==160)
H=np.linalg.inv(G)
dual=np.load(P/'dual_lll.npz')['U']
candidates=np.concatenate([np.eye(36,dtype=np.int64),dual])
norms=np.einsum('bi,ij,bj->b',candidates,H,candidates)
order=np.argsort(norms)
out=[]
for trial in order[:12]:
 c=candidates[trial]; norm=norms[trial]
 pivot=np.flatnonzero(np.abs(c)==1)[0]
 inds=[j for j in range(36) if j!=pivot]
 B=np.eye(36,dtype=np.int64)[inds]
 B[:,pivot]=-c[inds]*c[pivot]
 assert np.all(B@c==0)
 S=B@G@B.T; W=V[V@c==0][:,inds]
 # Work with one sign representative per pair.
 W=np.unique(np.where((W!=0).any(1)[:,None],W,0),axis=0)
 W=W[np.argmax(W!=0,axis=1).argsort(kind='stable')]
 signs=np.array([w[np.flatnonzero(w)[0]]>0 for w in W]);W=W[signs]
 # Whiten for stable LP scaling. In whitened coordinates the density derivative
 # is tr(E); x E x^T >= 0 preserves every known active length to first order.
 C=np.linalg.cholesky(S)
 X=W@C/np.sqrt(8)
 pairs=[(i,j) for i in range(35) for j in range(i,35)]
 A=np.column_stack([X[:,i]*X[:,j]*(1 if i==j else 2) for i,j in pairs])
 obj=np.array([float(i==j) for i,j in pairs])
 res=linprog(obj,A_ub=-A,b_ub=np.zeros(len(A)),bounds=[(-1,1)]*len(pairs),method='highs',options={'time_limit':90})
 r=dict(trial=int(trial),normal=c.tolist(),dual_norm=float(norm),inherited_contacts=int(2*len(W)),density_min8=float(16777216/3176523/np.sqrt(2*norm)),status=res.message,objective=float(res.fun) if res.fun is not None else None)
 print(json.dumps(r),flush=True);out.append(r)
 if res.success and res.fun<-1e-7:
  E=np.zeros((35,35))
  for (i,j),a in zip(pairs,res.x):E[i,j]=E[j,i]=a
  D=C@E@C.T
  np.savez(P/f'section{trial}_flex.npz',B=B,Q=20*S,D=D,W=W,normal=c,trace=res.fun)
  break
 (P/'contact_screen.json').write_text(json.dumps(out,indent=2))
(P/'contact_screen.json').write_text(json.dumps(out,indent=2))
