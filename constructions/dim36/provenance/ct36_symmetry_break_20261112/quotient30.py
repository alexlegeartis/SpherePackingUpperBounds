import numpy as np,json,itertools,time
from pathlib import Path
P=Path(__file__).resolve().parent;C=np.load(P.parent/'kp36_code/quaternary18.npz')['C'];w=np.count_nonzero(C,axis=1);rows=[]
def clique(adj):
 best=[]
 def go(cur,cand):
  nonlocal best
  if len(cur)+cand.bit_count()<=len(best):return
  if not cand:
   if len(cur)>len(best):best=cur
   return
  # greedy coloring upper bound
  vs=[];cols=[];left=cand;color=0
  while left:
   color+=1;avail=left
   while avail:
    b=avail&-avail;v=b.bit_length()-1;vs.append(v);cols.append(color);left^=b;avail^=b;avail&=~adj[v]
  for v,co in reversed(list(zip(vs,cols))):
   if len(cur)+co<=len(best):return
   go(cur+[v],cand&adj[v]);cand&=~(1<<v)
 go([0],adj[0]);return best
for deleted in [(0,1,2),(0,1,8),(8,9,17)]:
 pw=w-np.count_nonzero(C[:,deleted],axis=1)
 for extra in [x for x in range(18)if x not in deleted]:
  inds=deleted+(extra,);keys=np.sum(C[:,inds].astype(int)*(4**np.arange(4)),axis=1);m=np.full(256,100,dtype=int);np.minimum.at(m,keys[1:],pw[1:]);assert m[0]>=8
  allowed=m>=8;allowed[0]=False;adj=[]
  for i in range(256):adj.append(sum(1<<int(x)for x in np.flatnonzero(allowed[np.arange(256)^i])))
  cl=clique(adj);rows.append({'deleted':deleted,'extra':extra,'clique':cl,'size':len(cl)*1024})
 print(deleted,'max',max(r['size']for r in rows if r['deleted']==deleted),flush=True)
json.dump(rows,open(P/'quotient30.json','w'),indent=2)
