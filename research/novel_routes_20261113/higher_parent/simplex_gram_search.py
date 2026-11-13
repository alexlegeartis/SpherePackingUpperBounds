from pathlib import Path
import itertools,numpy as np,math,json
p=Path(__file__).resolve().parent
out=[]
for k in range(1,6):
 cand=[]
 if k==1:
  matrices=[np.array([[x]]) for x in range(8,34,2)]
 else:
  inds=list(itertools.combinations(range(k),2)); offvals=range(-4,1) if k<=4 else range(-3,1)
  # Fixed min8 diagonals; not a global optimization claim.
  matrices=(np.diag([8]*k)+sum((np.outer(np.eye(k,dtype=int)[i],np.eye(k,dtype=int)[j])+np.outer(np.eye(k,dtype=int)[j],np.eye(k,dtype=int)[i]))*v for (i,j),v in zip(inds,vs)) for vs in itertools.product(offvals,repeat=len(inds)))
 W=np.array(list(itertools.product(range(-1,2),repeat=k)),dtype=int);W=W[np.any(W,axis=1)]
 best=-1e99; count=0
 for K in matrices:
  count+=1
  sign,ld=np.linalg.slogdet(K)
  if sign<=0 or np.linalg.eigvalsh(K)[0]<1e-8:continue
  if np.einsum('bi,ij,bj->b',W,K,W).min()<8:continue
  H=np.linalg.inv(K);beta=max(float(np.diag(H).max()),float((np.diag(H)[:,None]+np.diag(H)[None,:]-2*H).max()));n=72-k
  if beta>=8:continue
  logd=math.log(k+1)+(n/2)*math.log((8-beta)/4)-ld/2
  if logd>best+1e-12:best=logd;bestK=K;bestbeta=beta
  cand.append((logd,K.tolist(),beta))
 cand.sort(key=lambda x:-x[0]);uniq=[];seen=set()
 for l,K,beta in cand:
  inv=np.linalg.inv(K); key=(round(l,10),round(beta,10),round(np.linalg.det(K)))
  if key in seen:continue
  seen.add(key);uniq.append({'K':K,'beta_approx':beta,'center_density_approx':math.exp(l),'log2_density':l/math.log(2),'det':round(np.linalg.det(K))})
  if len(uniq)>=20:break
 row={'rank':k,'dimension':72-k,'matrices_examined':count,'top':uniq};out.append(row)
 print(k,count,json.dumps(uniq[:3]),flush=True)
 (p/'simplex_gram_search.json').write_text(json.dumps(out,indent=2))
