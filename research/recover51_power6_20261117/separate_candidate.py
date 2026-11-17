"""Numerical all-face/edge/multistart separator, never a proof."""
import json,sys,itertools
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from scipy.special import eval_gegenbauer
from scipy.optimize import differential_evolution,minimize_scalar
D=Path(__file__).parent;path=Path(sys.argv[1]);v=json.load(open(path));n=v['dimension'];p=v['powers'][0];co=np.array(v['co']);B=np.array(v['B']);k=np.arange(len(co));den=eval_gegenbauer(k,(n-2)/2,1)
def value(x):
 x=np.asarray(x);a=1+6/13*x[...,0];b=a+6/13*(1-x[...,0])*x[...,1];c=b+6/13*(1-x[...,0])*(1-x[...,1])*x[...,2];u=(1-x[...,0])**4;w=u*(1-x[...,1])**4;z=w*(1-x[...,2])**4
 def bil(t,u,v):
  up=np.stack([np.ones_like(u),u,u*u],-1);vp=np.stack([np.ones_like(v),v,v*v],-1);G=eval_gegenbauer(k,(n-2)/2,t[...,None])/den*co
  return np.einsum('...k,kij,...i,...j->...',G,B,up,vp)
 return bil((a+b-c)/(2*np.sqrt(a*b)),u,w)+(1-x[...,2])**p*bil((a+c-b)/(2*np.sqrt(a*c)),u,z)+((1-x[...,1])*(1-x[...,2]))**p*bil((b+c-a)/(2*np.sqrt(b*c)),w,z)
pts=[];reports=[]
for axis in range(3):
 for fixed in [0.,1.]:
  def inject(y):return np.insert(y,axis,fixed)
  for seed in [17,51]:
   r=differential_evolution(lambda y:-value(inject(y)),[(0,1)]*2,seed=seed,maxiter=200,popsize=18,tol=1e-9);x=inject(r.x);pts.append(x);reports.append({'type':'face','axis':axis,'fixed':fixed,'point':x.tolist(),'value':float(value(x))})
 for fixed in itertools.product([0.,1.],repeat=2):
  grid=np.linspace(0,1,513);xx=np.tile(np.insert(fixed,axis,0.),(len(grid),1));xx[:,axis]=grid;vv=value(xx)
  for j in range(1,len(grid)-1):
   if vv[j]>=vv[j-1] and vv[j]>=vv[j+1]:
    def fun(t):x=xx[j].copy();x[axis]=t;return -value(x)
    r=minimize_scalar(fun,bounds=(grid[j-1],grid[j+1]),method='bounded',options={'xatol':1e-14});x=xx[j].copy();x[axis]=r.x;pts.append(x);reports.append({'type':'edge','point':x.tolist(),'value':float(value(x))})
for seed in [13,51,71]:
 r=differential_evolution(lambda x:-value(x),[(0,1)]*3,seed=seed,maxiter=250,popsize=20,tol=1e-9);pts.append(r.x);reports.append({'type':'interior','point':r.x.tolist(),'value':float(-r.fun)})
worst=max(reports,key=lambda r:r['value']);report={'status':'NUMERICAL_ONLY','candidate':path.name,'worst':worst,'reports':reports};(D/(path.stem+'_separation.json')).write_text(json.dumps(report,indent=2));cuts=D/f'guarded_cuts_p{p}.json';old=json.load(open(cuts))if cuts.exists()else[];cuts.write_text(json.dumps(old+[x.tolist()for x in pts]));print(worst,flush=True)
