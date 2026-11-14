"""Multiscale triangle separation. Only confirmed positive points are counterexamples."""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from scipy.optimize import minimize
from flint import fmpq as Q,arb,arb_poly,ctx
if not __debug__:raise RuntimeError('Assertions required')
ctx.prec=512
ap=argparse.ArgumentParser();ap.add_argument('kernel');ap.add_argument('--seeds',type=int,default=80)
args=ap.parse_args();file=Path(args.kernel).resolve();D=file.parent;r=json.loads(file.read_text())
n,d,m=r['dimension'],r['angular_degree'],r['features'];indices=list(map(tuple,r['angular_indices']))
rates=np.array([float(F(z)) for z in r['rates_exact']]);gamma=rates-rates[0]
bell=np.array([float(F(z)) for z in r['reference_coefficients']]);mat=np.array([[[float(F(z)) for z in row] for row in a] for a in r['blocks']])
coeff=np.stack([bell*mat[:,i,j] for i,j in indices],axis=1)
AP=[arb_poly([arb(Q(z)) for z in p]) for p in r['angular_polynomials']]
Tmax,Umax=float(F(r['compact_T_upper_exact'])),float(F(r['compact_u_upper_exact']))
def gegen(t):
    t=np.asarray(t);g=[np.ones_like(t),t]
    for ell in range(1,d):g.append(((2*ell+n-2)*t*g[-1]-ell*g[-2])/(ell+n-2))
    return np.stack(g,axis=-1)
def value(t,u,T):
    w2=1+u*u-2*t*u;w=np.sqrt(w2)
    rr=[np.exp(-(v[...,None]-1)*gamma) for v in [T,T*u*u,T*w2]]
    angles=[t,(1-u*t)/w,(u-t)/w]
    factors=[np.ones_like(t),np.exp(-rates[0]*T*(w2-u*u)),np.exp(-rates[0]*T*(w2-1))]
    out=np.zeros_like(t,dtype=float)
    for angle,(left,right),factor in zip(angles,[(0,1),(0,2),(1,2)],factors):
        angular=gegen(angle)@coeff
        for k,(i,j) in enumerate(indices):
            weight=rr[left][...,i]*rr[right][...,j]
            if i!=j:weight+=rr[left][...,j]*rr[right][...,i]
            out+=factor*weight*angular[...,k]
    return out
def exact(t,u,T):
    t,u,T=[arb(Q(str(F(float(z))))) for z in [t,u,T]]
    w2=1+u*u-2*t*u;w=w2.sqrt()
    rr=[[(-arb(Q(str(F(float(g)))))*(v-1)).exp() for g in gamma] for v in [T,T*u*u,T*w2]]
    ans=arb(0)
    for angle,(left,right),factor in zip([t,(1-u*t)/w,(u-t)/w],[(0,1),(0,2),(1,2)],
            [arb(1),(-arb(rates[0])*T*(w2-u*u)).exp(),(-arb(rates[0])*T*(w2-1)).exp()]):
        for k,(i,j) in enumerate(indices):
            weight=rr[left][i]*rr[right][j]
            if i!=j:weight+=rr[left][j]*rr[right][i]
            ans+=factor*weight*AP[k](angle)
    return ans
tt,uu,TT=np.meshgrid(np.r_[np.linspace(-1,0,81),np.linspace(0,.5,121)],
                  np.r_[np.linspace(1,1.3,51),np.linspace(1.3,2,21),np.linspace(2,Umax,11)],
                  np.r_[1,1.01,1.025,1.05,1.075,1.1,1.15,1.2,1.3,1.5,2,Tmax])
t,u,T=tt.ravel(),uu.ravel(),TT.ravel();mask=2*t*u<=1;t,u,T=t[mask],u[mask],T[mask]
vals=value(t,u,T);chosen=[];bins=set()
for j in np.argsort(vals)[::-1]:
    key=(int((t[j]+1)*100),int((u[j]-1)*50),int((T[j]-1)*100))
    if key in bins:continue
    bins.add(key);chosen.append(j)
    if len(chosen)>=args.seeds:break
counter=[];points=[];worst=-float('inf')
for j in chosen:
    def fun(z):
        tt,v,w=z;umax=min(Umax,1/(2*tt)) if tt>0 else Umax
        return -float(value(np.array(tt),np.array(1+v*(umax-1)),np.array(1+w*(Tmax-1))))
    umax=min(Umax,1/(2*t[j])) if t[j]>0 else Umax
    opt=minimize(fun,[t[j],(u[j]-1)/(umax-1) if umax>1 else 0,(T[j]-1)/(Tmax-1)],bounds=[(-1,.5),(0,1),(0,1)],method='L-BFGS-B',options={'ftol':1e-14,'gtol':1e-7,'maxiter':150})
    tt,v,w=opt.x;umax=min(Umax,1/(2*tt)) if tt>0 else Umax
    uu=1+v*(umax-1);TT=1+w*(Tmax-1)
    if tt>0 and 2*tt*uu>1:uu=np.nextafter(1/(2*tt),1.)
    val=exact(tt,uu,TT);worst=max(worst,float(val.mid()))
    point=[float(tt),float(uu),float(TT)];points.append({'point':point,'value_interval':str(val)})
    if val>0:counter.append(point)
cutsfile=D/'cuts.json';cuts=json.loads(cutsfile.read_text()) if cutsfile.exists() else {}
cuts['triangles']=cuts.get('triangles',[])+counter;cutsfile.write_text(json.dumps(cuts,indent=2)+'\n')
out={'certified_counterexamples':len(counter),'maximum_located_value':worst,'points':points,
     'scope':'Positive points verified rigorously; absence of positives is not a domain proof'}
(D/'numerical_separation.json').write_text(json.dumps(out,indent=2)+'\n')
print('MULTISCALE SEPARATION',len(counter),'maximum',worst,flush=True)
