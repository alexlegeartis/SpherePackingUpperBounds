from pathlib import Path
import numpy as np,json,itertools,time,math
from sympy.matrices.normalforms import smith_normal_form
import sympy as sp
p=Path(__file__).resolve().parent;G=np.load(p/'G72.npy');E=np.eye(72,dtype=np.int16)
V=[x for x in E]
for i in range(72):
 for j in range(i+1,72):
  if G[i,j]:
   v=E[i]-int(np.sign(G[i,j]))*E[j]
   if int(v@G@v)<=12:V.append(v)
V=np.unique(np.vstack((V,np.load(p/'gamma72_norm8_pool.npz')['V'])),axis=0);V=np.vstack((V,-V));VG=V.astype(np.int64)@G;norms=np.einsum('ij,ij->i',VG,V)
print('POOL',len(V),dict(zip(*np.unique(norms,return_counts=True))),flush=True)
K6=9*np.eye(6,dtype=int)-np.ones((6,6),dtype=int)
for i,j in [(0,1),(0,2),(1,3),(2,4),(3,5)]:K6[i,j]=K6[j,i]=-2
rows=[{'rank':6,'top':[{'K':K6.tolist()}]}];out=[]
cache={}
def ips(i):
 if i not in cache:
  cache[i]=VG@V[i]
  if len(cache)>128:cache.pop(next(iter(cache)))
 return cache[i]
def find(K):
 k=len(K);calls=0
 def rec(ids,options):
  nonlocal calls
  calls+=1
  if len(ids)==k:return ids
  lev=len(ids)
  for i in options:
   nxt=np.flatnonzero((norms==K[lev+1,lev+1]) & (ips(i)==K[lev,lev+1])) if lev+1<k else []
   if lev+1<k:
    for q,j in enumerate(ids):nxt=nxt[ips(j)[nxt]==K[q,lev+1]]
   if lev+1<k and len(nxt)==0:continue
   ans=rec(ids+[int(i)],nxt)
   if ans is not None:return ans
  return None
 ans=rec([],np.flatnonzero(norms==K[0,0]))
 return ans,calls
for row in rows:
 if row['rank']!=6:continue
 k=row['rank'];target=row['top'][0];K=np.array(target['K'],dtype=int);t=time.time();ids,calls=find(K);print('RESULT',k,ids,calls,time.time()-t,flush=True)
 if ids is not None:
  B=V[ids].astype(int);assert np.array_equal(B@G@B.T,K)
  snf=smith_normal_form(sp.Matrix(B),domain=sp.ZZ);diag=[abs(int(snf[i,i])) for i in range(k)];print('SNF',diag,flush=True)
  np.savez(p/f'gamma72_simplex_rank{k}.npz',G=G,B=B,K=K,S=np.vstack((np.zeros((1,k),dtype=int),np.eye(k,dtype=int))))
  out.append({'rank':k,'dimension':72-k,'normal_gram':K.tolist(),'normal_coefficients':B.tolist(),'saturation_diagonal':diag,**target})
(p/'embedded_simplex_rank6.json').write_text(json.dumps(out,indent=2))
