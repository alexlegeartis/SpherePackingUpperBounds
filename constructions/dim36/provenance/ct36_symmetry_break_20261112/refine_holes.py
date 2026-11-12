import numpy as np,json
from pathlib import Path
from fractions import Fraction as F
P=Path(__file__).resolve().parent;N=np.load(P/'reduced.npz')['N'];j=json.load(open(P/'holes.json'));V=np.loadtxt(P/'vectors.txt',dtype=np.int64).reshape(-1,36);VV=np.einsum('bi,ij,bj->b',V,N,V);assert np.all(VV==32)
out=[]
for r in j['unrejected']:
 p=r['p'];a=np.array(r['a']);y=a-p*np.array(r['witness']);k=0
 while k<30:
  g=V@(N@y);ix=int(np.argmax(abs(g)));sgn=1 if g[ix]>0 else -1;delta=p*p*32-2*p*abs(g[ix])
  if delta>=0:break
  y-=p*sgn*V[ix];k+=1
 q=F(int(y@N@y),4*p*p);ratio=F(24576,15625)*(min(q,F(8))/8)**18
 r.update({'refined_witness':((a-y)//p).tolist(),'refined_q':str(q),'refined_upper_ratio':str(ratio),'moves':k});out.append(r)
 print(p,q,'upperratio',float(ratio),flush=True)
json.dump(out,open(P/'refined_holes.json','w'),indent=2)
