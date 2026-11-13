"""Finite integer normal-Gram exploration for fixed nine-label clusters.
No candidate is certified unless embedded in Gamma72 and replayed separately.
"""
from pathlib import Path
import numpy as np,itertools,time,math,json
P=Path(__file__).resolve().parent;rng=np.random.default_rng(663399)
results=[]
for k in [6,7]:
 a=np.load(P/f'gamma72_root_cluster_rank{k}.npz');K0=a['K'];S=a['S'];W=np.array(list(itertools.product(range(-1,2),repeat=k)),dtype=np.int16);W=W[np.any(W,axis=1)];W=W[np.argmax(W!=0,axis=1)==0] if False else W
 W2=np.stack([np.outer(x,x).ravel() for x in W]);Z=np.array([x-y for i,x in enumerate(S) for y in S[:i]]);Z2=np.stack([np.outer(x,x).ravel() for x in Z]);n=72-k
 def score(K):
  q=W2@K.ravel()
  if q.min()<8:return -1e50
  sign,ld=np.linalg.slogdet(K)
  if sign<=0:return -1e50
  try:np.linalg.cholesky(K);H=np.linalg.inv(K)
  except np.linalg.LinAlgError:return -1e50
  beta=float((Z2@H.ravel()).max())
  if beta>=8:return -1e50
  return math.log(9)+n/2*math.log((8-beta)/4)-ld/2
 best=score(K0);K=K0.copy();old=best;bestK=K.copy();accepted=0;valid=0;start=time.time()
 for step in range(60000):
  KK=K.copy();i,j=map(int,rng.choice(k,2,replace=False));KK[i,j]+=int(rng.choice([-1,1]));KK[j,i]=KK[i,j]
  if step%5==0:
   ii,jj=map(int,rng.choice(k,2,replace=False));KK[ii,jj]+=int(rng.choice([-1,1]));KK[jj,ii]=KK[ii,jj]
  if np.max(np.abs(KK-np.diag(np.diag(KK))))>4:continue
  val=score(KK)
  if val<-1e20:continue
  valid+=1;temp=.03+.12*(1-(step%5000)/5000)
  if val>old or rng.random()<math.exp(min(0,(val-old)/temp)):K=KK;old=val;accepted+=1
  if val>best+1e-12:best=val;bestK=KK.copy();print('BEST',k,step,math.exp(best),bestK.tolist(),flush=True)
  if step%5000==4999:K=K0.copy();old=score(K)
 row={'rank':k,'attempted':60000,'screen_valid':valid,'accepted':accepted,'best_density_screen':math.exp(best),'bestK':bestK.tolist(),'starting_density':math.exp(score(K0)),'seconds':time.time()-start,'complete':False,'min_screen_only':True};results.append(row)
 print('DONE',row,flush=True);(P/'normal_mutation.json').write_text(json.dumps(results,indent=2))
