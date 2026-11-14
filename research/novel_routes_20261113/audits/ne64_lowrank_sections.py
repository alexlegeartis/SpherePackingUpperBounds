"""Search exact D5/E6/E7 normal embeddings; unrelated to floating SVP proof.

Returned integer embeddings inherit Ne64's published minimum12. Exact Gram,
saturation and density checks are performed for every saved construction.
"""
from pathlib import Path
import numpy as np,sympy as sp,json,time,itertools
from sympy.matrices.normalforms import smith_normal_form
P=Path(__file__).resolve().parent;BASE=P.parent/'dual64'
D=np.load(BASE/'reduced_dual.npz');G=D['G'];H=D['H'];T=D['T'];R=D['R']
V=np.load(BASE/'short_pool.npz')['V'].astype(np.int64)
patterns={'D5':(5,[(0,1),(1,2),(2,3),(2,4)]),'E6':(6,[(0,1),(1,2),(2,3),(3,4),(2,5)]),'E7':(7,[(0,1),(1,2),(2,3),(3,4),(4,5),(2,6)])}
rng=np.random.default_rng(461057);results=[]
def canonical(V):
 V=V.copy();V*=np.sign(V[np.arange(len(V)),np.argmax(V!=0,axis=1)])[:,None]
 return np.unique(V,axis=0)
def search(name,V,seconds=20):
 k,edges=patterns[name];K=12*np.eye(k,dtype=np.int64)
 for i,j in edges:K[i,j]=K[j,i]=-6
 degree=np.sum(K==-6,axis=1);order=sorted(range(k),key=lambda i:-degree[i]);Kp=K[np.ix_(order,order)]
 W=np.vstack([V,-V]);WH=W@R;cache={};start=time.time();calls=0
 def ips(j):
  if j not in cache:cache[j]=WH@W[j]
  return cache[j]
 def rec(ids,candidates):
  nonlocal calls
  calls+=1
  if calls%1000==0 and time.time()-start>seconds:return None
  d=len(ids)
  if d==k:return ids
  for j in candidates:
   if time.time()-start>seconds:return None
   if d+1<k:
    can=np.flatnonzero(ips(j)==Kp[d,d+1])
    for a,t in enumerate(ids):can=can[ips(t)[can]==Kp[a,d+1]]
    if not len(can):continue
   else:can=[]
   out=rec(ids+[int(j)],can)
   if out is not None:return out
  return None
 ids=rec([],rng.permutation(len(W)));print(name,'pool',len(V),'calls',calls,'time',time.time()-start,'found',ids,flush=True)
 if ids is None:return None
 Ap=W[ids]@T;A=Ap[np.argsort(order)];assert np.array_equal(A@H@A.T,K)
 S=smith_normal_form(sp.Matrix(A.tolist()),domain=sp.ZZ);diag=[abs(int(S[i,i]))for i in range(k)];assert diag==[1]*k
 det=sp.Matrix(K.tolist()).det()/3**k;n=64-k;delta_sq=sp.Rational(3**n,3**32)/det
 row={'name':name,'dimension':n,'normal_rank':k,'normal_det':str(det),'center_density_squared':str(delta_sq),'center_density':float(sp.sqrt(delta_sq)),'saturation_diagonal':diag,'ambient_minimum_input':12,'source_pool_size':len(V),'scope':'Exact construction. Global record status not asserted.'}
 np.savez_compressed(P/f'ne64_section_{name}.npz',G=G,H=H,A=A,Knum=K)
 return row
for level in range(3):
 assert np.all(np.einsum('bi,ij,bj->b',V,R,V)==12)
 for name in patterns:
  if any(r['name']==name for r in results):continue
  out=search(name,V,30)
  if out is not None:results.append(out)
 (P/'ne64_lowrank_sections.json').write_text(json.dumps(results,indent=2))
 if len(results)==3:break
 if level<2:
  gens=np.load(P.parent/'code_routes/cq64_generators.npz')['generators'][:5];Ti=np.array(sp.Matrix(T).inv(),dtype=np.int64);parts=[V]
  for g in gens:
   action=T@np.array(sp.Matrix(g).inv().T,dtype=np.int64)@Ti
   assert np.array_equal(action@R@action.T,R)
   parts.append(V@action)
  V=canonical(np.vstack(parts))
  if len(V)>20000:V=V[rng.choice(len(V),20000,False)]
  print('expanded',len(V),flush=True)
np.savez_compressed(P/'ne64_section_search_pool.npz',V=V.astype(np.int16))
print(json.dumps(results,indent=2),flush=True)
