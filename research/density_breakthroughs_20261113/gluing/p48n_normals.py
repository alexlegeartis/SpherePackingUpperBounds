"""Search rank9/10 primitive normal sublattices for39D/38D packing density.
All proposed records are checked by exact Gram determinant and Smith saturation.
The mother-lattice minimum6 is an explicitly inherited published input.
"""
from pathlib import Path
from collections import OrderedDict
import numpy as np,sympy as sp,json,math,time,argparse,re
from sympy.matrices.normalforms import smith_normal_form
P=Path(__file__).resolve().parent;ROOT=P.parents[2]
ap=argparse.ArgumentParser();ap.add_argument('--seconds',type=float,default=500);ap.add_argument('--seed',type=int,default=131310);ap.add_argument('--pool',default='P48n_glue_pool.npz');ap.add_argument('--tag',default='P48n');a=ap.parse_args()
j=np.load(P/a.pool);G=j['G'];V=j['V'].astype(np.int16);T=j['E8'].astype(np.int16);T*=np.sign(T[np.arange(8),np.argmax(T!=0,axis=1)])[:,None]
V=np.unique(V,axis=0);norms=np.einsum('bi,ij,bj->b',V.astype(np.int64),G,V.astype(np.int64));N=len(V);lookup={v.tobytes():i for i,v in enumerate(V)};seed8=[lookup[v.tobytes()]for v in T];VG=(V.astype(np.int64)@G).astype(np.float64);cache=OrderedDict()
def col(i):
 if i in cache:z=cache.pop(i);cache[i]=z;return z
 z=np.rint(VG@V[i].astype(float)).astype(np.int8);assert z[i]==norms[i];cache[i]=z
 if len(cache)>384:cache.popitem(last=False)
 return z
def score(ids):
 if not ids:return norms.astype(float).copy()
 C=np.column_stack([col(i)for i in ids]).astype(float);H=C[ids];inv=np.linalg.inv(H);q=norms-np.einsum('ij,ij->i',C@inv,C);q[q<1e-7]=np.inf;return q
start=time.time();rng=np.random.default_rng(a.seed);best={9:float('inf'),10:float('inf')};rawbest=best.copy();history=[];targets={9:21870,10:54675};seeds={};seen=set()
def save(ids,why):
 rank=len(ids);key=tuple(sorted(ids))
 if key in seen:return
 seen.add(key);B=V[ids].astype(np.int64);K=B@G@B.T;rawdet=int(sp.Matrix(K).det());assert rawdet>0
 S=smith_normal_form(sp.Matrix(B),domain=sp.ZZ);index=math.prod(abs(int(S[i,i]))for i in range(rank));assert index and rawdet%(index*index)==0;det=rawdet//(index*index)
 if det>=best[rank]:return
 best[rank]=det;dim=48-rank;density=1.5**(dim/2)/math.sqrt(det);r={'normal_rank':rank,'dimension':dim,'primitive_determinant':det,'raw_determinant':rawdet,'saturation_index':index,'center_density':density,'comparison_determinant':targets[rank],'density_ratio_to_comparison':math.sqrt(targets[rank]/det),'new_record_candidate':det<targets[rank],'pool_indices':list(map(int,ids)),'seconds':time.time()-start,'reason':why,'minimum_dependency':'published extremal P48n minimum6, not re-enumerated here'};history.append(r);np.savez(P/f'{a.tag}_best_rank{rank}.npz',G=G,B=B,K=K);(P/f'{a.tag}_normals_history.json').write_text(json.dumps(history,indent=2));print(json.dumps(r),flush=True)
print('POOL',N,'E8',seed8,flush=True)
for k in [9,10]:
 ids=seed8.copy()
 while len(ids)<k:
  q=score(ids);ids.append(int(rng.choice(np.flatnonzero(q<q.min()+1e-7))))
 seeds[k]=ids.copy();save(ids,'known_E8_extension_control')
trial=0
while time.time()-start<a.seconds:
 k=[9,10][trial%2]
 if trial%5==0:ids=seeds[k].copy()
 elif trial%5 in[1,2]:ids=[seed8[t]for t in sorted(rng.choice(8,int(rng.integers(3,7)),replace=False))]
 else:
  ids=[]
  for _ in range(int(rng.integers(1,5))):
   q=score(ids);ids.append(int(rng.choice(np.flatnonzero(np.isfinite(q)))))
 while len(ids)<k:
  q=score(ids);ids.append(int(rng.choice(np.flatnonzero(q<q.min()+1e-7))))
 for step in range(100):
  K=V[ids].astype(np.int64)@G@V[ids].astype(np.int64).T;d=float(np.linalg.det(K))
  if d<rawbest[k]-.01:rawbest[k]=d;save(ids,'determinant_improvement')
  elif step%12==0:save(ids,'primitive_saturation_check')
  count=int(rng.choice([1,1,2,2,3,4]));rem=set(map(int,rng.choice(k,count,replace=False)));ids=[x for i,x in enumerate(ids)if i not in rem]
  while len(ids)<k:
   q=score(ids);bestq=q.min();tied=np.flatnonzero(q<bestq+1e-7)
   if step%11==0 and len(tied)<50:tied=np.argsort(q)[:50]
   ids.append(int(rng.choice(tied)))
  if time.time()-start>=a.seconds:break
 trial+=1
 print('TRIAL',trial,'best',best,'seconds',time.time()-start,flush=True)
print('DONE',trial,best,time.time()-start,flush=True)
