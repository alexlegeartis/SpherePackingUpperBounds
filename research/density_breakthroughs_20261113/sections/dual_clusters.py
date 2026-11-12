"""Exact pair tests and bounded maximum-clique searches in a finite dual pool."""
from pathlib import Path
import numpy as np,time,json,sys
P=Path(__file__).resolve().parent
stem=sys.argv[1]if len(sys.argv)>1 else 'K11'
F=np.load(P/f'{stem}_extension.npz');A=F['adj'];det=int(F['det']);k=len(A);n=48-k
raw=np.loadtxt(P/f'{stem}_dual_vectors.txt',dtype=np.int64).reshape(-1,k+1)
X=np.concatenate([raw[:,1:],-raw[:,1:]])
qq=np.einsum('bi,ij,bj->b',X,A,X)
DD=qq[:,None]+qq[None,:]-2*(X@A@X.T)
bench={37:2**2.5,36:16777216/3176523,35:2**1.5,34:3**16.5/2**25}[n]
out=[]
def clique(adj,cap=20):
 best=[];start=time.time();calls=0;capped=False
 def colors(P):
  order=[];bounds=[];color=0;rem=P
  while rem:
   color+=1;Q=rem
   while Q:
    b=Q&-Q;v=b.bit_length()-1;order.append(v);bounds.append(color)
    rem^=b;Q^=b;Q&=~adj[v]
  return order,bounds
 def visit(C,P):
  nonlocal best,calls,capped
  calls+=1
  if calls%1024==0 and time.time()-start>cap:capped=True;return
  order,bounds=colors(P)
  for z in range(len(order)-1,-1,-1):
   if len(C)+bounds[z]<=len(best):return
   v=order[z];Q=P&adj[v]
   if Q:visit(C+[v],Q)
   elif len(C)+1>len(best):best=C+[v]
   P&=~(1<<v)
   if capped:return
 visit([], (1<<len(adj))-1)
 return best,calls,capped
for b in sorted(set(qq.tolist())):
 if b>det:continue
 idx=np.flatnonzero(qq<=b);D=DD[np.ix_(idx,idx)];ok=D<=b;np.fill_diagonal(ok,False)
 adj=[sum(1<<int(j)for j in np.flatnonzero(row))for row in ok]
 ids,calls,capped=clique(adj)
 Y=X[idx[ids]];m=len(Y)+1;beta=b/det
 delta=m*((6-beta)/4)**(n/2)/det**.5
 r={'rank':k,'dimension':n,'determinant':det,'beta_numerator':int(b),'beta_denominator':det,'size':m,'center_density':delta,'ratio':delta/bench,'vertices':len(idx),'clique_calls':calls,'capped':capped,'cluster':np.concatenate([np.zeros((1,k),dtype=np.int64),Y]).tolist()}
 assert np.all(D[np.ix_(ids,ids)]<=b)
 print(json.dumps({key:val for key,val in r.items()if key!='cluster'}),flush=True);out.append(r)
 (P/f'{stem}_clusters.json').write_text(json.dumps(out,indent=2))
