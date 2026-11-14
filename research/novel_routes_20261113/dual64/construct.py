"""Build antipode certificates using primitive sublattices of Ne64's dual.

All returned embeddings and cluster bounds are checked exactly. The ambient
minimum12 is inherited from the published catalogue lattice.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np,sympy as sp,json,itertools,time
from pathlib import Path
from sympy.matrices.normalforms import smith_normal_form
P=Path(__file__).resolve().parent
d=np.load(P/'reduced_dual.npz');G=d['G'];H=d['H'];T=d['T'];R=d['R'];n=len(G)
raw=np.loadtxt(P/'short_raw.txt',dtype=np.int64,ndmin=2)
V=raw[:,1:];assert np.all(np.einsum('bi,ij,bj->b',V,R,V)==12)
V=np.unique(V,axis=0).astype(np.int16);rng=np.random.default_rng(631361)
for step in range(8):
 ids=rng.choice(len(V),min(350,len(V)),False);C=V[ids].astype(np.int64)@R@V.T
 ii,jj=np.where(np.abs(C)==6);W=V[ids[ii]]-np.sign(C[ii,jj])[:,None]*V[jj]
 W=W.astype(np.int16);W*=np.sign(W[np.arange(len(W)),np.argmax(W!=0,axis=1)])[:,None]
 V=np.unique(np.vstack([V,W]),axis=0)
 if len(V)>40000:V=V[rng.choice(len(V),40000,False)]
 print('pool',step,len(V),flush=True)
 if len(V)==40000:break
assert np.all(np.einsum('bi,ij,bj->b',V.astype(np.int64),R,V.astype(np.int64))==12)
np.savez_compressed(P/'short_pool.npz',V=V)
V=np.vstack([V,-V]).astype(np.int64);VR=V@R;cache={}
def ips(i):
 if i not in cache:cache[i]=VR@V[i]
 return cache[i]
def embed(K):
 k=len(K);calls=0
 def rec(ids,options):
  nonlocal calls
  calls+=1
  if calls>100000:return None
  if len(ids)==k:return ids
  lev=len(ids)
  for j in options:
   if lev+1<k:
    nxt=np.flatnonzero(ips(j)==K[lev,lev+1])
    for u,i in enumerate(ids):nxt=nxt[ips(i)[nxt]==K[u,lev+1]]
    if not len(nxt):continue
   else:nxt=[]
   r=rec(ids+[int(j)],nxt)
   if r is not None:return r
  return None
 ids=rec([],np.arange(len(V)))
 return ids,calls
targets=[np.array([[12,-6],[-6,12]]),16*np.eye(3,dtype=int)-4*np.ones((3,3),dtype=int),15*np.eye(4,dtype=int)-3*np.ones((4,4),dtype=int)]
report=[]
def save(A,Knum,S,label):
 k=len(A);snf=smith_normal_form(sp.Matrix(A.tolist()),domain=sp.ZZ);diag=[abs(int(snf[i,i]))for i in range(k)]
 assert diag==[1]*k
 assert np.array_equal(A@H@A.T,Knum)
 K=sp.Matrix(Knum.tolist())/3;Ki=K.inv();beta=max(((sp.Matrix((x-y).tolist()).T*Ki*sp.Matrix((x-y).tolist()))[0]for x,y in itertools.combinations(S,2)),default=sp.Rational(0))
 dim=64-k;delta2=sp.Rational(len(S)**2,3**32)*((12-beta)/4)**dim/K.det()
 np.savez_compressed(P/f'{label}.npz',G=G,H=H,A=A,Knum=Knum,S=S)
 row={'label':label,'dimension':dim,'normal_rank':k,'cosets':len(S),'normal_gram_numerator':Knum.tolist(),'normal_gram_denominator':3,'normal_det':str(K.det()),'beta':str(beta),'delta_squared':str(delta2),'delta':float(sp.sqrt(delta2)),'saturation_diagonal':diag,'ambient_minimum_input':12,'status':'verified construction using published ambient minimum; global record status unestablished'}
 report.append(row);print(json.dumps(row),flush=True)
# The rank-one continuous optimum is dual norm16/3, available in this basis.
i=int(np.flatnonzero(R.diagonal()==16)[0]);save(T[[i]],np.array([[16]]),np.array([[0],[1]]),'ne64_rank1')
for K in targets:
 k=len(K);ids,calls=embed(K);print('embed',k,ids,calls,flush=True)
 if ids is not None:
  A=V[ids]@T;S=np.vstack([np.zeros(k,dtype=int),np.eye(k,dtype=int)])
  save(A,K,S,f'ne64_simplex_rank{k}')
(P/'constructed.json').write_text(json.dumps(report,indent=2)+'\n')
