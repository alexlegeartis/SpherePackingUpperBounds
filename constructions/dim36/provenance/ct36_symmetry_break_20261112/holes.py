import numpy as np,json
from pathlib import Path
from fractions import Fraction as F
P=Path(__file__).resolve().parent;N=np.load(P/'reduced.npz')['N'];R=np.linalg.cholesky(N).T; rng=np.random.default_rng(20261112)
rows=[];keep=[]
for p in (2,3,4):
 bestall=0
 for trial in range(300):
  a=rng.integers(0,p,36);best=None;bz=None
  for it in range(16):
   z=np.zeros(36,dtype=np.int64)
   for i in range(35,-1,-1):
    x=a[i]/p+R[i,i+1:]@(a[i+1:]/p-z[i+1:])/R[i,i]
    z[i]=round(x+(rng.normal(0,.23) if it else 0))
   y=a-p*z; grad=N@y
   for k in range(200):
    shifts=np.rint(-grad/(p*np.diag(N))).astype(np.int64);improve=2*p*shifts*grad+p*p*shifts*shifts*np.diag(N);i=int(np.argmin(improve))
    if improve[i]>=0:break
    y[i]+=p*shifts[i];grad+=p*shifts[i]*N[:,i]
   q=int(y@N@y)
   if best is None or q<best:best=q;bz=(a-y)//p
  f=F(best,4*p*p);ratio=F(24576,15625)*(f/8)**18
  rows.append({'p':p,'a':a.tolist(),'witness':bz.tolist(),'q':str(f),'upper_ratio':str(ratio)})
  bestall=max(bestall,float(f))
  if ratio>1:keep.append(rows[-1])
 print('p',p,'largest nearest witness',bestall,'notrejected',len(keep),flush=True)
json.dump({'rows':rows,'unrejected':keep,'scope':'900 sampled rational shifts; witnesses are exact upper bounds on hole distance'},open(P/'holes.json','w'))
