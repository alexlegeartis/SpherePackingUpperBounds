"""Sampled cyclic angular LP. Only candidate generation, never a sign proof."""
import numpy as np,json,math,argparse
from pathlib import Path
from scipy.optimize import linprog,minimize
D=Path(__file__).parent
ap=argparse.ArgumentParser();ap.add_argument('--dimension',type=int,default=45);ap.add_argument('--degree',type=int,default=19);args=ap.parse_args();n=args.dimension;deg=args.degree

def gegen(t):
 t=np.asarray(t);q=[np.ones_like(t),t]
 for l in range(1,deg):q.append(((2*l+n-2)*t*q[-1]-l*q[-2])/(l+n-2))
 return np.stack(q,axis=-1)
def cyc(t,u):
 w=np.sqrt(1+u*u-2*u*t)
 return gegen(t)+np.exp(-40*(1-2*u*t))[...,None]*gegen((1-u*t)/w)+np.exp(-40*(u*u-2*u*t))[...,None]*gegen((u-t)/w)
t=np.linspace(0,.5,251);u=np.r_[np.linspace(1,1.2,201),np.linspace(1.2,2,101),np.linspace(2,5,61)];T,U=np.meshgrid(t,u);T=T.ravel();U=U.ravel();mask=2*T*U<=1;T=T[mask];U=U[mask]
A=np.vstack([gegen(np.linspace(-1,.5001,6001)),cyc(T,U)]);b=np.full(len(A),-1e-4);rng=np.random.default_rng(20261111+n)
for rnd in range(10):
 res=linprog(np.ones(deg+1),A_ub=A,b_ub=b,A_eq=np.eye(deg+1)[0:1],b_eq=[1],bounds=[(0,None)]*(deg+1),method='highs',options={'dual_feasibility_tolerance':1e-8,'primal_feasibility_tolerance':1e-8})
 if not res.success:print(res.message);break
 c=res.x;ext=[]
 for t0 in np.linspace(-1,.5001,501):
  lo=max(-1,t0-.004);hi=min(.5001,t0+.004)
  z=minimize(lambda v:-float(gegen(v[0])@c),[t0],bounds=[(lo,hi)],method='L-BFGS-B')
  if -z.fun>-1e-5:ext.append(gegen(z.x[0]))
 vals=cyc(T,U)@c;ind=np.argsort(vals)[-60:]
 for j in ind:
  # Parameter v in[0,1] imposes u∈[1,min(5,1/(2t))].
  def fun(z):
   tt,v=z;umax=min(5,1/(2*tt)) if tt>0 else 5;uu=1+v*(umax-1)
   return -float(cyc(np.array(tt),np.array(uu))@c)
  umax=min(5,1/(2*T[j])) if T[j]>0 else 5;v0=(U[j]-1)/(umax-1) if umax>1 else 0
  z=minimize(fun,[T[j],v0],bounds=[(0,.5),(0,1)],method='L-BFGS-B')
  if -z.fun>-1e-5:
   tt,v=z.x;umax=min(5,1/(2*tt)) if tt>0 else 5;ext.append(cyc(np.array(tt),np.array(1+v*(umax-1))))
 out={'dimension':n,'degree':deg,'coefficients':c.tolist(),'B':float(res.fun),'round':rnd,'extra':len(ext),'certified':False}
 print(out,flush=True);(D/f'angular{n}_candidate.json').write_text(json.dumps(out,indent=2))
 if not ext:break
 A=np.vstack([A,np.array(ext)]);b=np.r_[b,np.full(len(ext),-1e-3)]
