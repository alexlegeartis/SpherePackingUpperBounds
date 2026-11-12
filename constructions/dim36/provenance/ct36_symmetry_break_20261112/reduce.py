import numpy as np
from pathlib import Path
P=Path(__file__).resolve().parent
G=np.load(P/'N.npy');T=np.eye(36,dtype=np.int64);H=G.copy()
def gs(H):
 n=len(H);mu=np.zeros((n,n));d=np.zeros(n)
 for i in range(n):
  for j in range(i):mu[i,j]=(H[i,j]-sum(mu[i,k]*mu[j,k]*d[k]for k in range(j)))/d[j]
  d[i]=H[i,i]-sum(mu[i,k]**2*d[k]for k in range(i))
 return mu,d
k=1;s=0
while k<36:
 mu,d=gs(H)
 for j in range(k-1,-1,-1):
  q=round(mu[k,j])
  if q:T[k]-=q*T[j];H=T@G@T.T;mu,d=gs(H)
 if d[k]>=(.99-mu[k,k-1]**2)*d[k-1]:k+=1
 else:T[[k,k-1]]=T[[k-1,k]];H=T@G@T.T;k=max(1,k-1)
 s+=1
np.savez(P/'reduced.npz',N=H,T=T);np.savetxt(P/'R.txt',np.linalg.cholesky(H/4).T);print('LLL',s)
