"""Numerical 3D separation followed by rigorous checks of proposed counterexamples."""
import argparse,json,sys
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from scipy.optimize import minimize
from flint import fmpq as Q,arb,arb_poly,ctx
from common import rational
if not __debug__:raise RuntimeError('Assertions required')
ctx.prec=512
ap=argparse.ArgumentParser();ap.add_argument('kernel');ap.add_argument('--seeds',type=int,default=60)
args=ap.parse_args();file=Path(args.kernel).resolve();D=file.parent;r=json.loads(file.read_text())
n=r['dimension'];d=r['angular_degree'];slow=float(F(r['slow_exact']));gamma=float(F(r['fast_exact']))-slow
Umax=float(F(r['compact_u_upper_exact']));Tmax=float(F(r['compact_T_upper_exact']))
bell=np.array([float(F(z)) for z in r['reference_coefficients']]);mats=np.array([[[float(F(z)) for z in row] for row in m] for m in r['blocks']])
AP=[arb_poly([arb(Q(c)) for c in p]) for p in r['angular_polynomials']]
def gegen(t):
    t=np.asarray(t);g=[np.ones_like(t),t]
    for k in range(1,d):g.append(((2*k+n-2)*t*g[-1]-k*g[-2])/(k+n-2))
    return np.stack(g,axis=-1)
def angular(t):
    g=gegen(t)
    return tuple(g@(bell*mats[:,i,j]) for i,j in [(0,0),(0,1),(1,1)])
def value(t,u,T):
    w2=1+u*u-2*t*u;w=np.sqrt(w2)
    rA=np.exp(-gamma*(T-1));rB=np.exp(-gamma*(T*u*u-1));rC=np.exp(-gamma*(T*w2-1))
    def term(z,x,y):
        p,q,r=angular(z);return p+(x+y)*q+x*y*r
    return term(t,rA,rB)+np.exp(-slow*T*(w2-u*u))*term((1-u*t)/w,rA,rC)+np.exp(-slow*T*(w2-1))*term((u-t)/w,rB,rC)
def exact(t,u,T):
    t,u,T=[arb(Q(str(F(float(z))))) for z in [t,u,T]]
    w2=1+u*u-2*t*u;w=w2.sqrt()
    rA=(-arb(gamma)*(T-1)).exp();rB=(-arb(gamma)*(T*u*u-1)).exp();rC=(-arb(gamma)*(T*w2-1)).exp()
    def term(z,x,y):return AP[0](z)+(x+y)*AP[1](z)+x*y*AP[2](z)
    return term(t,rA,rB)+(-arb(slow)*T*(w2-u*u)).exp()*term((1-u*t)/w,rA,rC)+(-arb(slow)*T*(w2-1)).exp()*term((u-t)/w,rB,rC)
t,u,T=np.meshgrid(np.r_[np.linspace(-1,0,101),np.linspace(0,.5,151)],
                  np.r_[np.linspace(1,1.3,61),np.linspace(1.3,2,31),np.linspace(2,Umax,21)],
                  np.r_[1.,1.01,1.025,1.05,1.075,1.1,1.15,1.2,1.3,1.5,2.,Tmax])
t,u,T=t.ravel(),u.ravel(),T.ravel();mask=2*t*u<=1;t,u,T=t[mask],u[mask],T[mask]
vals=value(t,u,T);indices=np.argsort(vals)[-args.seeds:]
counter=[];worst=-float('inf');points=[]
for j in indices:
    def fun(z):
        tt,v,ww=z;umax=min(Umax,1/(2*tt)) if tt>0 else Umax
        return -float(value(np.array(tt),np.array(1+v*(umax-1)),np.array(1+ww*(Tmax-1))))
    umax=min(Umax,1/(2*t[j])) if t[j]>0 else Umax
    x0=[t[j],(u[j]-1)/(umax-1) if umax>1 else 0,(T[j]-1)/(Tmax-1)]
    opt=minimize(fun,x0,bounds=[(-1,.5),(0,1),(0,1)],method='L-BFGS-B',options={'ftol':1e-14,'gtol':1e-7,'maxiter':150})
    tt,v,ww=opt.x;umax=min(Umax,1/(2*tt)) if tt>0 else Umax
    uu=1+v*(umax-1);TT=1+ww*(Tmax-1)
    if tt>0 and 2*tt*uu>1:uu=np.nextafter(1/(2*tt),1.)
    val=exact(tt,uu,TT);worst=max(worst,float(val.mid()))
    points.append({'point':[float(tt),float(uu),float(TT)],'value_interval':str(val)})
    if val>0:counter.append([float(tt),float(uu),float(TT)])
cutsfile=D/'cuts.json';cuts=json.loads(cutsfile.read_text()) if cutsfile.exists() else {}
cuts['triangles']=cuts.get('triangles',[])+counter
cutsfile.write_text(json.dumps(cuts,indent=2)+'\n')
out={'certified_counterexamples':len(counter),'maximum_located_value':worst,
     'points':points,'scope':'Counterexamples rigorous; absence of counterexamples is not a full-domain proof'}
(D/'numerical_separation.json').write_text(json.dumps(out,indent=2)+'\n')
print('MATRIX SEPARATION',len(counter),'maximum located',worst,flush=True)
