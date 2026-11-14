"""Section exchange using all 167096 sampled minima, with cached dot columns."""
from pathlib import Path
from collections import OrderedDict
import numpy as np,sympy as s,re,math,json,time,sys
from sympy.matrices.normalforms import smith_normal_form
P=Path(__file__).resolve().parent;ROOT=P.parents[2]
seconds=float(sys.argv[1])if len(sys.argv)>1 else 300
F=np.load(P/'P48p_pool.npz');G=F['G'];V=np.concatenate([F['V'],np.load(P/'P48p_enumeration/pool.npz')['V']]).astype(np.int16)
sign=np.sign(V[np.arange(len(V)),np.argmax(V!=0,axis=1)]);V*=sign[:,None]
V=np.unique(V,axis=0);N=len(V);lookup={v.tobytes():j for j,v in enumerate(V)}
VG=(V.astype(np.int64)@G).astype(float)
cert=ROOT/'research/density_breakthroughs_20261113/seed_deformations/external_sphere_packing/packings_out/records/dim38-delta9/P48p-a17f250f.cert.txt'
text=cert.read_text();a=text.index('[[',text.index('coefficient rows'));b=text.index(']]',a)+2
T=np.array([[int(z)for z in re.findall('-?\\d+',l)]for l in text[a:b].splitlines()],dtype=np.int16)
sg=np.sign(T[np.arange(len(T)),np.argmax(T!=0,axis=1)]);T*=sg[:,None]
seed=[lookup[v.tobytes()]for v in T]
cache=OrderedDict()
def col(i):
 if i in cache:
  z=cache.pop(i);cache[i]=z;return z
 z=np.rint(VG@V[i].astype(float)).astype(np.int8);cache[i]=z
 if len(cache)>384:cache.popitem(last=False)
 return z
def scores(ids):
 if not ids:return np.full(N,6.)
 C=np.column_stack([col(i)for i in ids]).astype(float)
 H=C[ids];inv=np.linalg.inv(H)
 q=6-np.einsum('ij,ij->i',C@inv,C);q[q<1e-7]=np.inf
 return q
start=time.time();rng=np.random.default_rng(13700);best={11:121500,12:301725,13:692334,14:1412199};rawbest=best.copy();history=[]
bench={11:2**2.5,12:16777216/3176523,13:2**1.5,14:3**16.5/2**25}
def save(ids):
 k=len(ids);B=V[ids].astype(np.int64);H=B@G@B.T;det=int(s.Matrix(H).det())
 S=smith_normal_form(s.Matrix(B),domain=s.ZZ);idx=math.prod(abs(int(S[i,i]))for i in range(k));assert det%(idx*idx)==0;det//=idx*idx
 if det>=best[k]:return
 best[k]=det;n=48-k;delta=1.5**(n/2)/math.sqrt(det)
 row={'rank':k,'dimension':n,'determinant':det,'saturation_index':idx,'center_density':delta,'ratio':delta/bench[k],'elapsed':time.time()-start}
 print(json.dumps(row),flush=True);history.append(row);np.savez(P/f'large_best_k{k}.npz',G=G,B=B,K=H)
 (P/'large_section_history.json').write_text(json.dumps(history,indent=2))
trial=0
print('POOL',N,flush=True)
while time.time()-start<seconds:
 k=[11,13,14,12][trial%4]
 if trial%3==0:ids=seed.copy()
 elif trial%3==1:
  # Delete a small but nontrivial portion of E8/K10 before completion.
  ids=[seed[i]for i in sorted(rng.choice(10,5,replace=False))]
 else:ids=[int(rng.integers(N))]
 while len(ids)<k:
  q=scores(ids);x=q.min();ids.append(int(rng.choice(np.flatnonzero(q<x+1e-6))))
 for step in range(60):
  H=V[ids].astype(np.int64)@G@V[ids].astype(np.int64).T;d=float(np.linalg.det(H))
  if d<rawbest[k]-.01:rawbest[k]=d;save(ids)
  elif step%15==0:save(ids)
  remove=int(rng.choice([1,1,1,2,3]));r=set(map(int,rng.choice(len(ids),remove,False)));ids=[x for i,x in enumerate(ids)if i not in r]
  while len(ids)<k:
   q=scores(ids);x=q.min()
   cand=np.flatnonzero(q<x+1e-6)
   if step%7==0:cand=np.argsort(q)[:min(24,N)]
   ids.append(int(rng.choice(cand)))
  if time.time()-start>=seconds:break
 trial+=1
print('DONE','trials',trial,'best',best,'elapsed',time.time()-start,flush=True)
