"""Exchange higher-rank normal subspaces using exact P48n norm6/8/10/12 vectors.

The explicit ambient pool is produced by the gluing research directory.  All
reported determinants are exact and divided by the squared saturation index.
Minimum6 of the ambient lattice remains the published construction input.
"""
from pathlib import Path
from collections import OrderedDict
import numpy as np, sympy as s, math, json, time, sys
from sympy.matrices.normalforms import smith_normal_form
P=Path(__file__).resolve().parent
seconds=float(sys.argv[1])if len(sys.argv)>1 else 180
F=np.load(P.parent/'gluing/P48n_extended_normal_pool.npz')
G=F['G'];V=F['V'];norms=F['norms'];T=F['E8'].astype(np.int16)
sg=np.sign(T[np.arange(len(T)),np.argmax(T!=0,axis=1)]);T*=sg[:,None]
seed=[]
for v in T:
 ix=np.flatnonzero(np.all(V==v,axis=1));assert len(ix)==1;seed.append(int(ix[0]))
N=len(V); VG=(V.astype(np.int64)@G).astype(float)
assert np.all(np.einsum('ij,ij->i',VG,V)==norms)
cache=OrderedDict()
def col(i):
 if i in cache:
  z=cache.pop(i);cache[i]=z;return z
 z=np.rint(VG@V[i].astype(float)).astype(np.int8);cache[i]=z
 if len(cache)>256:cache.popitem(last=False)
 return z
def scores(ids):
 if not ids:return norms.astype(float)
 C=np.column_stack([col(i)for i in ids]).astype(float)
 H=C[ids];inv=np.linalg.inv(H)
 q=norms-np.einsum('ij,ij->i',C@inv,C);q[q<1e-7]=np.inf
 return q
start=time.time();rng=np.random.default_rng(91382)
best={11:10**20,13:10**20,14:10**20};history=[]
bench={11:2**2.5,13:2**1.5,14:3**16.5/2**25}
def save(ids):
 k=len(ids);B=V[ids].astype(np.int64);H=B@G@B.T
 S=smith_normal_form(s.Matrix(B),domain=s.ZZ)
 idx=math.prod(abs(int(S[i,i]))for i in range(k))
 rawdet=round(np.linalg.det(H));candidate=rawdet/(idx*idx)
 if candidate>=best[k]+.01:return
 det=int(s.Matrix(H).det());assert det%(idx*idx)==0;det//=idx*idx
 if det>=best[k]:return
 best[k]=det;n=48-k;delta=1.5**(n/2)/math.sqrt(det)
 row={'rank':k,'dimension':n,'determinant':det,'saturation_index':idx,'center_density':delta,'ratio':delta/bench[k],'elapsed':time.time()-start,'basis_norms':np.diag(H).tolist()}
 print(json.dumps(row),flush=True);history.append(row)
 np.savez(P/f'P48n_extended_best_k{k}.npz',G=G,B=B,K=H)
 (P/'P48n_extended_section_history.json').write_text(json.dumps(history,indent=2))
trial=0
print('POOL',N,flush=True)
while time.time()-start<seconds:
 k=[11,13,14][trial%3]
 if trial%4==0:ids=seed.copy()
 elif trial%4 in (1,2):ids=[seed[i]for i in sorted(rng.choice(8,rng.choice([3,5,6]),replace=False))]
 else:ids=[int(rng.integers(N))]
 while len(ids)<k:
  q=scores(ids);x=q.min();ids.append(int(rng.choice(np.flatnonzero(q<x+1e-6))))
 for step in range(35):
  save(ids)
  remove=int(rng.choice([1,1,1,2,3,4]));r=set(map(int,rng.choice(len(ids),remove,False)));ids=[x for i,x in enumerate(ids)if i not in r]
  while len(ids)<k:
   q=scores(ids);x=q.min();cand=np.flatnonzero(q<x+1e-6)
   if step%8==0:cand=np.argpartition(q,24)[:24]
   ids.append(int(rng.choice(cand)))
  if time.time()-start>=seconds:break
 trial+=1
print('DONE','trials',trial,'best',best,'elapsed',time.time()-start,flush=True)
