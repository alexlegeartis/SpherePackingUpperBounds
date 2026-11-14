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
V=np.unique(V,axis=0);V=np.vstack((V,-V));VG=V.astype(np.int64)@G;norms=np.einsum('ij,ij->i',VG,V)
print('POOL',len(V),dict(zip(*np.unique(norms,return_counts=True))),flush=True)
rows=json.loads((p/'simplex_gram_search.json').read_text());out=[]
cache={}
def ips(i):
 if i not in cache:cache[i]=VG@V[i]
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
 k=row['rank'];target=row['top'][0];K=np.array(target['K'],dtype=int);t=time.time();ids,calls=find(K);print('RESULT',k,ids,calls,time.time()-t,flush=True)
 if ids is not None:
  B=V[ids].astype(int);assert np.array_equal(B@G@B.T,K)
  snf=smith_normal_form(sp.Matrix(B),domain=sp.ZZ);diag=[abs(int(snf[i,i])) for i in range(k)];print('SNF',diag,flush=True)
  np.savez(p/f'gamma72_simplex_rank{k}.npz',G=G,B=B,K=K,S=np.vstack((np.zeros((1,k),dtype=int),np.eye(k,dtype=int))))
  out.append({'rank':k,'dimension':72-k,'normal_gram':K.tolist(),'normal_coefficients':B.tolist(),'saturation_diagonal':diag,**target})
(p/'embedded_simplex.json').write_text(json.dumps(out,indent=2))
