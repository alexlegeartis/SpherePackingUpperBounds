"""Finite heuristic dual anticodes for 4D5,4E6,4E7,4E8 normals in Gamma72.
Floating enumeration is discovery only. Every saved finite cluster is checked exactly.
"""
from pathlib import Path
import numpy as np,sympy as sp,json,math,time
P=Path(__file__).resolve().parent;rng=np.random.default_rng(723264)
def shortvectors(H,bound):
 k=len(H);R=np.linalg.cholesky(H).T;z=np.zeros(k,dtype=np.int64);out=[]
 def rec(i,rem):
  if i<0:out.append(z.copy());return
  center=-R[i,i+1:]@z[i+1:]/R[i,i];rad=math.sqrt(max(rem,0))/R[i,i];lo=math.ceil(center-rad-1e-8);hi=math.floor(center+rad+1e-8)
  for a in range(lo,hi+1):
   z[i]=a;rec(i-1,rem-(R[i,i]*(a-center))**2)
 rec(k-1,bound+1e-7);return np.array(out,dtype=np.int64)
results=[]
for k in range(5,9):
 f=P/f'gamma72_root_rank{k}.npz';a=np.load(f);K=a['K'];H=sp.Matrix(K.tolist()).inv();det=int(sp.Matrix(K.tolist()).det());J=np.array((H*det).tolist(),dtype=np.int64);W=shortvectors(np.array(H.tolist(),dtype=float),1.01);q=np.einsum('bi,ij,bj->b',W,J,W);W=W[q<=det];q=q[q<=det];origin=int(np.flatnonzero(np.all(W==0,axis=1))[0]);A=W@J@W.T;D=q[:,None]+q[None,:]-2*A;assert np.min(D)==0
 print('POOL',k,len(W),flush=True);best=2**((72-k)/2)/math.sqrt(det);bestS=[origin];bestP=0
 for beta in np.unique(D):
  if beta==0 or beta>det:continue
  inds=np.flatnonzero(q<=beta);DD=D[np.ix_(inds,inds)];adj=DD<=beta;zero=int(np.flatnonzero(inds==origin)[0]);degrees=adj.sum(axis=1);seeds=np.argsort(-degrees)[:min(40,len(inds))]
  for seed in seeds:
   chosen=[zero];opts=np.flatnonzero(adj[zero]);opts=opts[opts!=zero]
   if seed!=zero:chosen.append(int(seed));opts=opts[adj[seed,opts] & (opts!=seed)]
   while len(opts):
    subdeg=adj[np.ix_(opts,opts)].sum(axis=1);v=int(opts[np.argmax(subdeg)]);chosen.append(v);opts=opts[adj[v,opts] & (opts!=v)]
   S=inds[chosen];diam=int(D[np.ix_(S,S)].max());val=len(S)*((8-diam/det)/4)**((72-k)/2)/math.sqrt(det)
   if val>best+1e-5:best=val;bestS=S.copy();bestP=diam;print('BEST',k,len(S),diam/det,val,flush=True)
 data={'normal_rank':k,'dimension':72-k,'det':det,'cosets':len(bestS),'beta':str(sp.Rational(bestP,det)),'density':best,'Ne64_record_when_64':3**16 if k==8 else None,'beats_Ne64':best>3**16 if k==8 else None,'heuristic_discovery':True,'global_cluster_optimality':False}
 np.savez(P/f'gamma72_root_cluster_rank{k}.npz',G=a['G'],B=a['B'],K=K,S=W[bestS]);results.append(data);(P/'root_cluster_search.json').write_text(json.dumps(results,indent=2))
 print('DONE',data,flush=True)
