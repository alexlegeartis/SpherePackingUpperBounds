import numpy as np,sympy as sp,json
from pathlib import Path
P=Path(__file__).resolve().parent;B=np.load(P/'split_basis.npy');I=2*sp.Matrix(B).inv();assert all(x.q==1 for x in I);I=np.array(I,dtype=np.int64);piv={};inds=[]
for j,row in enumerate(I%2):
 z=sum(int(x)<<k for k,x in enumerate(row));w=z
 while w:
  k=w.bit_length()-1
  if k in piv:w^=piv[k]
  else:piv[k]=w;inds.append(j);break
assert len(inds)==12
C=np.zeros((4096,36),dtype=np.int64)
for k,i in enumerate(inds):C[(np.arange(4096)>>k)&1==1]^=I[i]%2
assert len(np.unique(C,axis=0))==4096
D=np.load(P/'flex.npz');G4=B@D['G4']@B.T;F=B@D['D']@B.T;assert not G4[:12,12:].any();assert not F[:12,12:].any();assert not F[:12,:12].any()
np.savez(P/'glue4096.npz',B=B,glue_twice=C,G4=G4,D=F,generator_indices=np.array(inds));print('PASS4096cosets blockorthogonality',flush=True)
