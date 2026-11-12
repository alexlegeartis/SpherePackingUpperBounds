import numpy as np,sympy as sp,math
D=np.load('/tmp/ct36/data.npz');H=D['H'];V=np.load('/tmp/ct36/min6.npy');A=V.T@V@H
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp
S,U,W=smith_normal_decomp(DomainMatrix.from_Matrix(sp.Matrix(A)).convert_to(sp.ZZ))
from sympy.matrices.normalforms import hermite_normal_form
B=np.array(hermite_normal_form(W.to_Matrix()[:,24:]).T.tolist(),dtype=np.int64);
G=B@H@B.T;T=np.eye(12,dtype=np.int64)
exec(open('/tmp/ct36/prepare.py').read().split('def gs(H):')[1].split('k=1;')[0].join(['def gs(H):','']))
k=1
while k<12:
 F=T@G@T.T;mu,d=gs(F)
 for j in range(k-1,-1,-1):
  q=round(mu[k,j])
  if q:T[k]-=q*T[j];F=T@G@T.T;mu,d=gs(F)
 if d[k]>=(.99-mu[k,k-1]**2)*d[k-1]:k+=1
 else:T[[k,k-1]]=T[[k-1,k]];k=max(k-1,1)
print((T@G@T.T).diagonal());np.save('/tmp/ct36/kernel_sat_basis.npy',T@B)
