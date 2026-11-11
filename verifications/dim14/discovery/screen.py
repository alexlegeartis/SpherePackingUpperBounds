"""Numerical discovery only: curved separable radial damping, dimension 14.
All outcomes are diagnostics until exact PSD/Fourier/sign/cyclic proofs pass.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import argparse,json,math
from pathlib import Path
import numpy as np
from scipy.optimize import linprog,minimize
from scipy.special import eval_genlaguerre,gammaln
from numpy.polynomial import Polynomial
D=Path(__file__).resolve().parent
N=14; L=19

def gegen(t):
 t=np.asarray(t);q=[np.ones_like(t),t]
 for l in range(1,L):q.append(((2*l+N-2)*t*q[-1]-l*q[-2])/(l+N-2))
 return np.stack(q,axis=-1)

def cyc(t,u,a,b):
 w2=1+u*u-2*u*t;w=np.sqrt(w2)
 y=u*u;z=w2
 r1=np.exp(-a*(z-y)-b*((z-1)**2-(y-1)**2))
 r2=np.exp(-a*(z-1)-b*(z-1)**2)
 return gegen(t)+r1[...,None]*gegen((1-u*t)/w)+r2[...,None]*gegen((u-t)/w)

def angular(a,b,rounds=2):
 t=np.linspace(-1,.5,451);u=np.r_[np.linspace(1,1.3,121),np.linspace(1.3,2,51),np.linspace(2,5,31)]
 T,U=np.meshgrid(t,u);T=T.ravel();U=U.ravel();m=2*T*U<=1;T=T[m];U=U[m]
 A=np.vstack([gegen(np.linspace(-1,.5,5001)),cyc(T,U,a,b)]);rhs=np.full(len(A),-1e-6)
 ref=json.loads((D/'inputs/scaling_reference15.json').read_text())
 ss=np.maximum(np.asarray(ref['coefficients']),1.);ss[0]=1.
 G=[Polynomial([1]),Polynomial([0,1])];X=G[1]
 for l in range(1,L):G.append(((2*l+N-2)*X*G[-1]-l*G[-2])/(l+N-2))
 history=[]
 for rnd in range(rounds):
  r=linprog(ss/max(ss),A_ub=A*ss[None,:],b_ub=rhs,A_eq=np.eye(L+1)[:1],b_eq=[1],bounds=[(0,None)]*(L+1),method='highs',options={'time_limit':90,'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
  if not r.success:return {'success':False,'reason':r.message,'a':a,'b':b}
  c=r.x*ss;P=sum((v*g for v,g in zip(c,G)),Polynomial([0]));roots=P.deriv().roots();ts=[-1.,.5]+[z.real for z in roots if abs(z.imag)<1e-7 and -1<z.real<.5]
  ext=[gegen(t0) for t0 in ts if P(t0)>-1e-6]
  vals=cyc(T,U,a,b)@c
  # Numerical boundary and all near-top interior separation, not a proof.
  starts=np.r_[np.argsort(vals)[-40:], np.where(U==1)[0][np.argsort(vals[U==1])[-10:]]]
  worst=float(vals.max())
  for j in starts:
   def fn(z):
    t0,v=z;umax=min(5,1/(2*t0)) if t0>0 else 5;uu=1+v*(umax-1)
    return -float(cyc(np.asarray(t0),np.asarray(uu),a,b)@c)
   umax=min(5,1/(2*T[j])) if T[j]>0 else 5;v=(U[j]-1)/(umax-1) if umax>1 else 0
   opt=minimize(fn,[T[j],v],method='L-BFGS-B',bounds=[(-1,.5),(0,1)],options={'maxiter':100,'ftol':1e-14,'gtol':1e-10})
   worst=max(worst,float(-opt.fun))
   if -opt.fun>-1e-7:
    t0,v=opt.x;umax=min(5,1/(2*t0)) if t0>0 else 5;ext.append(cyc(np.asarray(t0),np.asarray(1+v*(umax-1)),a,b))
  # Root diagnostic for the proposed all-scale lemma; never accept numeric roots as proof.
  rr=P.roots();real=sorted(float(z.real) for z in rr if abs(z.imag)<1e-6 and .5<float(z.real)<1)
  history.append({'round':rnd,'B':float(sum(c)),'added':len(ext),'worst_sampled':worst,'roots_above_half':real})
  if not ext:break
  A=np.vstack([A,ext]);rhs=np.r_[rhs,np.full(len(ext),-1e-5)]
 return {'success':True,'dimension':N,'angular_degree':L,'a':a,'b':b,'B':float(sum(c)),'coefficients':c.tolist(),'history':history,'certified':False}

def radial(a,b,B,degree=63,s=7/4,cap_profile=None):
 d=degree;cut=40.
 norms=np.exp(gammaln(np.arange(d+1)+N/2)-gammaln(N/2)-gammaln(np.arange(d+1)+1));sign=(-1.)**np.arange(d+1)
 def lag(X):return np.array([eval_genlaguerre(k,N/2-1,2*X)/norms[k] for k in range(d+1)]).T
 y=np.r_[np.linspace(0,4*d+N+20,6001),np.geomspace(4*d+N+20,2500,801)]
 x=np.unique(np.r_[np.linspace(np.pi*s,cut,1501),np.pi*s+np.geomspace(1e-8,2,501)])
 z=np.r_[np.linspace(cut,4*d+N+20,4001),np.geomspace(4*d+N+20,2500,801)]
 v=x/(np.pi*s)-1
 h=np.exp(-a*v-b*v*v) if cap_profile is None else cap_profile(v+1)
 g=np.exp(x)*(2*h-h*h)/B
 A=np.vstack([np.c_[lag(x),-g],np.c_[lag(z),np.zeros(len(z))],np.c_[-lag(y)*sign,np.zeros(len(y))]])
 scale=np.maximum(abs(A).max(axis=1),1e-200);A/=scale[:,None]
 eq=np.r_[sign,0.];c=np.ones(d+2)
 r=linprog(c,A_ub=A,b_ub=np.full(len(A),-1e-9),A_eq=eq[None,:],b_eq=[1.],bounds=[(None,None)]*(d+1)+[(0,None)],method='highs',options={'time_limit':90,'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
 if not r.success:return {'success':False,'reason':r.message}
 return {'success':True,'density':float(r.fun*np.pi**(N/2)*(s/4)**(N/2)/math.gamma(N/2+1)),'objective':float(r.fun),'C':float(r.x[-1]),'radial':r.x[:-1].tolist(),'max_scaled_residual':float(max(A@r.x)),'normalization_residual':float(eq@r.x-1),'degree':d,'s':s,'cut_X':cut,'certified':False}

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--a',type=float,required=True);ap.add_argument('--b',type=float,required=True);ap.add_argument('--rounds',type=int,default=2);args=ap.parse_args()
 assert args.a>0 and args.b>=0
 name=f'a{args.a:g}_b{args.b:g}';p=D/(name+'.json')
 out=angular(args.a,args.b,args.rounds);p.write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items() if k!='coefficients'},flush=True)
 if out['success']:
  out['radial']=radial(args.a,args.b,out['B']);p.write_text(json.dumps(out,indent=2)+'\n');print('RADIAL',name,{k:v for k,v in out['radial'].items() if k!='radial'},flush=True)
