"""Sampled cyclic angular LP. Only candidate generation, never a sign proof."""
import numpy as np,json,math,argparse
from pathlib import Path
from scipy.optimize import linprog,minimize
ROOT=Path(__file__).resolve().parent;D=ROOT/'work';D.mkdir(exist_ok=True)
from fractions import Fraction as F
from flint import fmpq_poly,fmpq,arb,ctx
ctx.prec=192
X=fmpq_poly([0,1])
ap=argparse.ArgumentParser();ap.add_argument('--dimension',type=int,default=45);ap.add_argument('--degree',type=int,default=19);ap.add_argument('--alpha',type=int,default=40);ap.add_argument('--rounds',type=int,default=3);args=ap.parse_args();n=args.dimension;deg=args.degree;alpha=args.alpha

def gegen(t):
 t=np.asarray(t);q=[np.ones_like(t),t]
 for l in range(1,deg):q.append(((2*l+n-2)*t*q[-1]-l*q[-2])/(l+n-2))
 return np.stack(q,axis=-1)
def cyc(t,u):
 w=np.sqrt(1+u*u-2*u*t)
 return gegen(t)+np.exp(-alpha*(1-2*u*t))[...,None]*gegen((1-u*t)/w)+np.exp(-alpha*(u*u-2*u*t))[...,None]*gegen((u-t)/w)
t=np.linspace(0,.5,251);u=np.r_[np.linspace(1,1.2,201),np.linspace(1.2,2,101),np.linspace(2,5,61)];T,U=np.meshgrid(t,u);T=T.ravel();U=U.ravel();mask=2*T*U<=1;T=T[mask];U=U[mask]
A=np.vstack([gegen(np.linspace(-1,.5001,6001)),cyc(T,U)]);b=np.full(len(A),-1e-4);rng=np.random.default_rng(20261111+n)
ref=json.loads((ROOT/'inputs/angular_scaling_reference.json').read_text());ss=np.array([float(F(x)) for x in ref['gegenbauer_coefficients'][:deg+1]]);ss[ss==0]=1e8;ss[0]=1.;ss=np.maximum(ss,1.)
for rnd in range(args.rounds):
 res=linprog(ss/max(ss),A_ub=A*ss[None,:],b_ub=b,A_eq=np.eye(deg+1)[0:1],b_eq=[1],bounds=[(0,None)]*(deg+1),method='highs',options={'dual_feasibility_tolerance':1e-8,'primal_feasibility_tolerance':1e-8})
 if not res.success:print(res.message);break
 c=res.x*ss;ext=[]
 G=[fmpq_poly([1]),X]
 for ell in range(1,deg):G.append(((2*ell+n-2)*X*G[-1]-ell*G[-2])/(ell+n-2))
 P=sum((fmpq_poly([fmpq(F(float(z)).numerator,F(float(z)).denominator)])*g for z,g in zip(c,G)),fmpq_poly())
 critical=[-1.,.5001]+[float(z.real.mid()) for z,m in P.derivative().numer().complex_roots() if z.imag.contains(0) and z.real>-1 and z.real<arb('0.5001')]
 for tt in critical:
  val=float(P(fmpq(F(tt).numerator,F(tt).denominator)))
  if val>-1e-4:ext.append(gegen(tt))
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
 out={'alpha':alpha,'dimension':n,'degree':deg,'coefficients':c.tolist(),'B':float(res.fun*max(ss)),'round':rnd,'extra':len(ext),'certified':False}
 print(out,flush=True);(D/f'angular{n}_a{alpha}_scaled_candidate.json').write_text(json.dumps(out,indent=2))
 if not ext:break
 A=np.vstack([A,np.array(ext)]);b=np.r_[b,np.full(len(ext),-1e-3)]
