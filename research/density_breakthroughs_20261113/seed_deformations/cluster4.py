"""Bounded finite search of dual anticodes, not a global minimum/embedding proof."""
import numpy as np,itertools,json,pathlib,time,sys
HERE=pathlib.Path(__file__).resolve().parent
record=509.6189826428642;best=record
T=np.array(list(itertools.product(range(-2,3),repeat=4)),dtype=float);T=T[np.any(T,axis=1)]
small=np.array(list(itertools.product(range(-1,2),repeat=4)),dtype=float);small=small[np.any(small,axis=1)]
ii,jj=np.triu_indices(4,1)

def clique(A,k):
 # Exact combinatorial branching on the supplied floating-screen graph.
 N=len(A);adj=[]
 for row in A:
  adj.append(sum(1<<int(i) for i in np.flatnonzero(row)))
 def go(cand,path):
  if len(path)>=k:return path
  if cand.bit_count()<k-len(path):return None
  while cand:
   b=cand&-cand;cand-=b;v=b.bit_length()-1
   z=go(cand&adj[v],path+[v])
   if z:return z
   if cand.bit_count()<k-len(path):return None
  return None
 return go((1<<N)-1,[])

hits=[];count=0;graphs=0;tic=time.time()
# Sign/permutation gauge reduction is only applied for equal diagonal6.
for diag in itertools.combinations_with_replacement([6,8,10,12],4):
 if diag<(6,6,12,12):continue
 # Gauge first row nonpositive; all remaining offdiagonals in [-4,4].
 for abc in itertools.product(range(-4,1),repeat=3):
  if diag==(6,6,6,6) and abc!=tuple(sorted(abc)):continue
  for tail in itertools.product(range(-4,5),repeat=3):
   off=abc+tail;G=np.diag(np.array(diag,dtype=float));G[ii,jj]=off;G[jj,ii]=off
   count+=1
   # Short explicit lattice vectors reject too-small minimum early.
   if np.min(np.einsum('ij,jk,ik->i',small,G,small))<5.999999:continue
   det=np.linalg.det(G)
   if det<323.999999 or np.linalg.eigvalsh(G)[0]<=.05:continue
   F=np.linalg.inv(G);qn=np.einsum('ij,jk,ik->i',T,F,T)
   for s in range(6,17):
    beta=6-4*(best*np.sqrt(det)/s)**(1/22)
    if beta<=0:continue
    S=T[qn<beta-1e-10]
    if len(S)<s-1:continue
    qs=np.einsum('ij,jk,ik->i',S,F,S);P=S@F@S.T
    adj=(qs[:,None]+qs[None,:]-2*P<beta-1e-10);np.fill_diagonal(adj,False)
    if np.count_nonzero(np.sum(adj,axis=1)>=s-2)<s-1:continue
    graphs+=1;ids=clique(adj,s-1)
    if ids:
     C=np.vstack([np.zeros(4),S[ids]]);ds=np.einsum('ij,jk,ik->i',C,F,C);actual=float((ds[:,None]+ds[None,:]-2*C@F@C.T).max());delta=s*((6-actual)/4)**22/np.sqrt(det)
     r={'s':s,'beta':actual,'delta':float(delta),'det':float(det),'G':G.astype(int).tolist(),'cluster':C.astype(int).tolist(),'minimum_screen_only':True,'embedding_verified':False}
     span=np.ceil(np.sqrt(np.maximum(0,6*np.diag(F)))).astype(int); vv=np.array(list(itertools.product(*[range(-int(a),int(a)+1) for a in span]))); vv=vv[np.any(vv,axis=1)]; qq=np.einsum('ij,jk,ik->i',vv,G,vv)
     if qq.min()<5.999999:continue
     hits.append(r);best=delta;print('HIT',json.dumps(r),flush=True);(HERE/'cluster4_hits.json').write_text(json.dumps(hits,indent=2))
 print('diag',diag,'count',count,'graphs',graphs,'best',best,'elapsed',time.time()-tic,flush=True)
 (HERE/'cluster4_progress.json').write_text(json.dumps({'count':count,'graphs':graphs,'best':best,'diag':diag,'elapsed':time.time()-tic,'hits':hits},indent=2))
