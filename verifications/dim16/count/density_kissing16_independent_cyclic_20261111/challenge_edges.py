if not __debug__:
    raise RuntimeError("Certificate verification requires assertions; do not run Python with -O or -OO.")

from challenge import *
from scipy.optimize import minimize_scalar
start=time.monotonic();results=[]
for s in [-1,1]:
 u=np.linspace(-1,C,1000001);v=u*C+s*np.sqrt(np.maximum(0,(1-u*u)*(1-C*C)));valid=(v>=-1)&(v<=C);u=u[valid];v=v[valid]
 vals=cycle(np.c_[u,v,np.full_like(u,C)])
 ids=np.r_[0,len(u)-1,np.where((vals[1:-1]>=vals[:-2])&(vals[1:-1]>=vals[2:]))[0]+1];ids=sorted(ids,key=lambda i:-vals[i])[:100]
 def cv(x):return float(cycle([x,x*C+s*np.sqrt(max(0,(1-x*x)*(1-C*C))),C]))
 best=(-np.inf,None)
 for i in ids:
  opt=minimize_scalar(lambda x:-cv(x),bounds=(u[max(0,i-10)],u[min(len(u)-1,i+10)]),method='bounded',options={'xatol':1e-15})
  if -opt.fun>best[0]:best=(-opt.fun,opt.x)
 val,x=best;pt=[x,x*C+s*np.sqrt(max(0,(1-x*x)*(1-C*C))),C]
 rec={'sign':s,'sampled':len(u),'max':val,'point':pt,'seconds':time.monotonic()-start};results.append(rec);print(rec,flush=True)
(ROOT/'challenge_edges.json').write_text(json.dumps({'scope':'Numerical Gram-boundary/cap-edge challenge; not proof.','results':results},indent=2))
