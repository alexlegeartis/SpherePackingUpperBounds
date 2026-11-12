"""Finite-pool random greedy section search; exact records checked on save."""
import os,time,json,sys,math
from pathlib import Path
import numpy as np
import sympy as s
from sympy.matrices.normalforms import smith_normal_form
P=Path(__file__).resolve().parent
ambient=sys.argv[1]if len(sys.argv)>1 else 'P48p'
seconds=float(sys.argv[2])if len(sys.argv)>2 else 300
data=np.load(P/f'{ambient}_pool.npz');V=data['V'];G=data['G'];N=len(V)
inn=np.rint((V@G).astype(float)@V.T).astype(np.int8)
assert np.all(np.diag(inn)==6)
rng=np.random.default_rng(2026111300);start=time.time()
bench={11:2**2.5,12:16777216/3176523,13:2**1.5,14:3**16.5/2**25}
best={k:float('inf')for k in bench}; rawbest=best.copy();history=[];seen=set()
def save(ids,reason):
 k=len(ids);B=V[ids];H=B@G@B.T
 det=int(s.Matrix(H).det())
 if det<=0:return
 S=smith_normal_form(s.Matrix(B),domain=s.ZZ)
 content=math.prod(abs(int(S[i,i]))for i in range(k))
 assert content and det%(content*content)==0
 det//=content*content
 if det>=best[k]:return
 best[k]=det;n=48-k;delta=1.5**(n/2)/math.sqrt(det)
 r={'ambient':ambient,'normal_rank':k,'dimension':n,'determinant':det,'unsaturated_determinant':int(s.Matrix(H).det()),'saturation_index':content,'center_density':delta,'benchmark':bench[k],'ratio':delta/bench[k],'reason':reason,'pool_indices':list(map(int,ids)),'elapsed':time.time()-start}
 print(json.dumps(r),flush=True);history.append(r)
 np.savez(P/f'{ambient}_best_k{k}.npz',G=G,B=B,gram=H)
 (P/f'{ambient}_section_history.json').write_text(json.dumps(history,indent=2))
def scores(ids):
 if not ids:return np.full(N,6.)
 H=inn[np.ix_(ids,ids)].astype(float)
 inv=np.linalg.inv(H)
 C=inn[:,ids].astype(float)
 q=6-np.einsum('ij,ij->i',C@inv,C)
 q[q<1e-7]=float('inf')
 return q
trial=0
while time.time()-start<seconds:
 k=[11,13,14,12][trial%4];ids=[]
 # Seed from varied four-vectors before greedy completion, so E8 is not forced.
 nseed=1 if trial%5==0 else rng.integers(2,6)
 for _ in range(nseed):
  q=scores(ids);valid=np.flatnonzero(np.isfinite(q));ids.append(int(rng.choice(valid)))
 for _ in range(k-len(ids)):
  q=scores(ids);bestq=q.min();tied=np.flatnonzero(q<bestq+1e-6)
  ids.append(int(rng.choice(tied)))
 for step in range(80):
  H=inn[np.ix_(ids,ids)].astype(float);det=float(np.linalg.det(H))
  if det<rawbest[k]-.01:
   rawbest[k]=det;save(ids,'raw_determinant_improvement')
  elif step%20==0:
   key=(k,tuple(sorted(ids)))
   if key not in seen:
    seen.add(key);save(ids,'saturation_screen')
  count=1 if step%7 else int(rng.choice([2,3]))
  removed=rng.choice(len(ids),count,replace=False)
  ids=[v for j,v in enumerate(ids)if j not in removed]
  for _ in range(count):
   q=scores(ids);bestq=q.min()
   if step%11==0:
    # Occasional warmer jump; no score here is accepted as a proof.
    order=np.argsort(q)[:32];idx=int(rng.choice(order))
   else:
    tied=np.flatnonzero(q<bestq+1e-6);idx=int(rng.choice(tied))
   ids.append(idx)
  if time.time()-start>=seconds:break
 trial+=1
print('DONE',ambient,'trials',trial,'best',best,flush=True)
