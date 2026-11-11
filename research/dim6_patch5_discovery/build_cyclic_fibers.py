"""Build numerical discovery rows from exact rational cyclic boundary polynomials.
These rows are not themselves a density proof; full exact replay is mandatory.
"""
import math,time,numpy as np
from pathlib import Path
from flint import fmpq as Q,fmpq_poly as P,fmpq_mat as M
X=P([0,1]);c=1+X/2;K=12;p=4;N=p+1;D=3*K;n=6;t0=time.time();T=[P([1]),2*X-1]
for j in range(1,p):T.append(2*T[1]*T[-1]-T[-2])
G=[P([1]),X]
for k in range(1,K):G.append(((2*k+n-2)*X*G[-1]-k*G[-2])/Q(k+n-2))
trans=[]
for h in range(16):
 lo=Q(h,16);w=Q(1,16)
 trans.append(M([[sum((Q(math.comb(j,m)*math.comb(i,m),math.comb(D,m))*lo**(j-m)*w**m for m in range(min(i,j)+1)),Q(0))for j in range(D+1)]for i in range(D+1)]))
rows=[]
for a,b in [(Q(1),Q(1)),(Q(1),Q(3,2)),(Q(3,2),Q(3,2))]:
 ta=[q(2*(a-1))for q in T];tb=[q(2*(b-1))for q in T]
 allp=[]
 for k in range(K+1):
  F=G[k]((a*a+b*b-c*c)/(2*a*b))*c**K
  A=sum((G[k][m]/(2*a)**m*(a*a+c*c-b*b)**m*c**(K-m)for m in range(k+1)),P([]))
  B=sum((G[k][m]/(2*b)**m*(b*b+c*c-a*a)**m*c**(K-m)for m in range(k+1)),P([]))
  for i in range(N):
   for j in range(N):allp.append(ta[i]*tb[j]*F+ta[i]*T[j]*A+tb[i]*T[j]*B)
 C=M([[q[j]for q in allp]for j in range(D+1)])
 for BT in trans:
  out=BT*C;rows.extend([[float(out[i,j])*math.exp(-math.pi)for j in range(len(allp))]for i in range(D+1)])
 print('FIBER',str(a),str(b),'seconds',time.time()-t0,flush=True)
arr=np.array(rows);dest=Path(__file__).with_name('cyclic_fiber_rows.npz');np.savez_compressed(dest,rows=arr);print('SAVED',arr.shape,'seconds',time.time()-t0)
