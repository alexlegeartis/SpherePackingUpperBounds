import numpy as np,re,sympy as sp,time,json
s=open('/tmp/ct36/source.html').read();s=s.split('NAME="GRAM"')[1].split('</pre>')[0];a=list(map(int,re.findall(r'-?\d+',s)));print(a[:5],len(a))
G=np.array(a[2:1298],dtype=np.int64).reshape(36,36);assert np.array_equal(G,G.T)
print('det',sp.Matrix(G).det());assert sp.Matrix(G).det()==3**18
T=np.eye(36,dtype=np.int64);H=G.copy()
def gs(H):
 n=len(H);mu=np.zeros((n,n));d=np.zeros(n)
 for i in range(n):
  for j in range(i):mu[i,j]=(H[i,j]-sum(mu[i,k]*mu[j,k]*d[k] for k in range(j)))/d[j]
  d[i]=H[i,i]-sum(mu[i,k]**2*d[k] for k in range(i))
 return mu,d
k=1;steps=0
while k<36:
 mu,d=gs(H)
 for j in range(k-1,-1,-1):
  q=round(mu[k,j])
  if q:T[k]-=q*T[j];H=T@G@T.T;mu,d=gs(H)
 if d[k]>=(.99-mu[k,k-1]**2)*d[k-1]:k+=1
 else:T[[k,k-1]]=T[[k-1,k]];H=T@G@T.T;k=max(k-1,1)
 steps+=1
 if steps>10000:raise Exception('lll limit')
print('lll',steps,'diagonal',H.diagonal());np.savez('/tmp/ct36/data.npz',G=G,H=H,T=T);np.savetxt('/tmp/ct36/R.txt',np.linalg.cholesky(H).T);json.dump({'det':str(3**18),'lll_steps':steps},open('/tmp/ct36/setup.json','w'))
