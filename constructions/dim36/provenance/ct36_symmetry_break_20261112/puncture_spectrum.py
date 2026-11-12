import numpy as np,json,itertools
from pathlib import Path
from fractions import Fraction as F
P=Path(__file__).resolve().parent;C=np.load(P.parent/'kp36_code/quaternary18.npz')['C'];N=len(C);w=np.count_nonzero(C,axis=1);out=[];rng=np.random.default_rng(42)
for p in (1,2,3):
 pools=list(itertools.combinations(range(18),p));sel=pools if p==1 else [pools[k]for k in rng.choice(len(pools),12,replace=False)]
 for deleted in sel:
  pw=w-np.count_nonzero(C[:,deleted],axis=1);edge=(pw<8)&(pw>0);a=edge.astype(np.int64);deg=int(a.sum())
  for s in [1<<i for i in range(18)]:
   b=a.reshape(-1,2*s);u=b[:,:s].copy();v=b[:,s:].copy();b[:,:s]=u+v;b[:,s:]=u-v
  lo=int(a.min());upper=F(-N*lo,deg-lo)
  out.append({'deleted':deleted,'length':18-p,'degree':deg,'eigen_min':lo,'independence_upper':str(upper)})
 print('puncture',p,'upper range',min(float(F(r['independence_upper']))for r in out if r['length']==18-p),max(float(F(r['independence_upper']))for r in out if r['length']==18-p),flush=True)
json.dump(out,open(P/'puncture_spectrum.json','w'),indent=2)
