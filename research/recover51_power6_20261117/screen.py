import json, numpy as np
from scipy.special import eval_gegenbauer
from scipy.optimize import differential_evolution
from fractions import Fraction
from pathlib import Path
D=Path(__file__).parent
v=json.load(open(D/'base/geometry/dim51/frozen.json'))
co=np.array([float(Fraction(x)) for x in v['co']]); B=np.array([[[float(Fraction(x))for x in row]for row in b]for b in v['B']]); n=v['dimension']; k=np.arange(len(co));den=eval_gegenbauer(k,(n-2)/2,1)
def value(x,p):
 x=np.asarray(x); a=1+6/13*x[...,0];b=a+6/13*(1-x[...,0])*x[...,1];c=b+6/13*(1-x[...,0])*(1-x[...,1])*x[...,2]
 u=(1-x[...,0])**4;w=u*(1-x[...,1])**4;z=w*(1-x[...,2])**4
 def bil(t,u,v):
  up=np.stack([np.ones_like(u),u,u*u],-1);vp=np.stack([np.ones_like(v),v,v*v],-1)
  G=eval_gegenbauer(k,(n-2)/2,t[...,None])/den*co
  return np.einsum('...k,kij,...i,...j->...',G,B,up,vp)
 return bil((a+b-c)/(2*np.sqrt(a*b)),u,w)+(1-x[...,2])**p*bil((a+c-b)/(2*np.sqrt(a*c)),u,z)+((1-x[...,1])*(1-x[...,2]))**p*bil((b+c-a)/(2*np.sqrt(b*c)),w,z)
rng=np.random.default_rng(51); xs=np.r_[rng.random((20000,3)),np.array(np.meshgrid(*[[0,.001,.01,.1,.5,.9,.99,.999,1]]*3)).reshape(3,-1).T];res=[]
for p in [4,5,6,7,8,9,10,12]:
 y=value(xs,p);ix=np.argmax(y);de=differential_evolution(lambda x:-value(x,p),[(0,1)]*3,seed=51,maxiter=120,popsize=12,tol=1e-9)
 rec={'power':p,'max_random':float(y[ix]),'random_point':xs[ix].tolist(),'max_optimized':float(-de.fun),'optimized_point':de.x.tolist(),'status':'NUMERICAL_ONLY'};res.append(rec);print(rec,flush=True);(D/'screen_results.json').write_text(json.dumps(res,indent=2))
