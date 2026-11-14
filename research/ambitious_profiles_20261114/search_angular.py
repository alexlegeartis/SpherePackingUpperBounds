"""Sampled angular LP for h(v)=v^-a exp[-b(v-1)]; numerical discovery only."""
import argparse,json,math,time
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog,minimize
from flint import fmpq as Q, fmpq_poly as P, arb,ctx
from common import HERE,base_angular,normalized_gegenbauer,rational

ap=argparse.ArgumentParser()
ap.add_argument('dimension',type=int);ap.add_argument('--degree',type=int,default=19)
ap.add_argument('--a',default='20');ap.add_argument('--b',default='20')
ap.add_argument('--rounds',type=int,default=5);ap.add_argument('--tail-u',type=int,default=8)
ap.add_argument('--tag',default='',help='Separate output directory for a retained follow-up search')
args=ap.parse_args();n,d=args.dimension,args.degree
a,b=F(args.a),F(args.b);assert a+b>0 and b>0
ctx.prec=256
tag=args.tag or f'n{n}_a{str(a).replace("/","_")}_b{str(b).replace("/","_")}_ang{d}'
D=HERE/tag;D.mkdir(exist_ok=True)
run=len(list(D.glob('round*')))
def gegen(t):
    t=np.asarray(t);g=[np.ones_like(t),t]
    for l in range(1,d):g.append(((2*l+n-2)*t*g[-1]-l*g[-2])/(l+n-2))
    return np.stack(g[:d+1],axis=-1)
def cyc(t,u):
    w2=1+u*u-2*u*t;w=np.sqrt(w2)
    r1=np.exp(-float(a)*np.log(w2/(u*u))-float(b)*(w2-u*u))
    r2=np.exp(-float(a)*np.log(w2)-float(b)*(w2-1))
    return gegen(t)+r1[...,None]*gegen((1-u*t)/w)+r2[...,None]*gegen((u-t)/w)
t=np.r_[np.linspace(-1,0,161),np.linspace(0,.5,181)]
u=np.r_[np.linspace(1,1.2,121),np.linspace(1.2,2,81),np.linspace(2,args.tail_u,81)]
T,U=np.meshgrid(t,u);T=T.ravel();U=U.ravel();mask=2*T*U<=1;T=T[mask];U=U[mask]
A=np.vstack([gegen(np.linspace(-1,.5001,8001)),cyc(T,U)])
rhs=np.full(len(A),-1e-4)
ref=json.loads(base_angular(n).read_text());sc=np.array([float(F(v)) for v in ref['gegenbauer_coefficients']])
sc=np.pad(sc,(0,max(0,d+1-len(sc))),constant_values=max(sc))[:d+1]
sc=np.maximum(sc,1.);sc[0]=1
G=normalized_gegenbauer(n,d);start=time.time()
for rnd in range(args.rounds):
    res=linprog(sc/max(sc),A_ub=A*sc[None,:],b_ub=rhs,A_eq=np.eye(d+1)[0:1],b_eq=[1],
                bounds=[(0,None)]*(d+1),method='highs',
                options={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9,'time_limit':180})
    if not res.success:raise RuntimeError(res.message)
    c=res.x*sc;c[0]=1
    cs=[Q(str(F(float(v)))) for v in c]
    eps=Q(1,1000);cs[0]-=eps;cs=[v/(1-eps) for v in cs]
    poly=sum((z*g for z,g in zip(cs,G)),P([]))
    exact={'dimension':n,'degree':d,'profile_a_exact':str(a),'profile_b_exact':str(b),
           'cutoff':'5001/10000','negative_cap_delta':'1/100000','tail_u':args.tail_u,
           'gegenbauer_coefficients':list(map(str,cs)),'polynomial_coefficients':list(map(str,poly.coeffs())),
           'B_exact':str(sum(cs,Q(0))),'constant_guard_exact':'1/1000','certified':False}
    numeric={'dimension':n,'degree':d,'a_exact':str(a),'b_exact':str(b),'certified':False,
             'raw_B':float(sum(c)),'guarded_B':float(sum(cs,Q(0))),'coefficients':c.tolist(),
             'round':rnd,'run':run,'seconds':time.time()-start}
    folder=D/f'round{run:03}_{rnd:02}';folder.mkdir()
    (folder/'candidate.json').write_text(json.dumps(numeric,indent=2)+'\n')
    (folder/'angular.json').write_text(json.dumps(exact,indent=2)+'\n')
    (D/'angular.json').write_text(json.dumps(exact,indent=2)+'\n')
    (D/'candidate.json').write_text(json.dumps(numeric,indent=2)+'\n')
    extra=[]
    rawpoly=sum((Q(str(F(float(z))))*g for z,g in zip(c,G)),P([]))
    for z,m in rawpoly.derivative().numer().complex_roots():
        if z.imag.contains(0) and z.real>arb(-1) and z.real<arb(Q(5001,10000)):
            tt=float(z.real.mid())
            if rawpoly(Q(str(F(tt))))>-Q(1,10000):extra.append(gegen(tt))
    vals=cyc(T,U)@c
    for j in np.argsort(vals)[-40:]:
        def fun(z):
            tt,v=z;umax=min(args.tail_u,1/(2*tt)) if tt>0 else args.tail_u
            uu=1+v*(umax-1)
            return -float(cyc(np.array(tt),np.array(uu))@c)
        umax=min(args.tail_u,1/(2*T[j])) if T[j]>0 else args.tail_u
        v0=(U[j]-1)/(umax-1) if umax>1 else 0
        z=minimize(fun,[T[j],v0],bounds=[(-1,.5),(0,1)],method='L-BFGS-B')
        if -z.fun>-1e-4:
            tt,v=z.x;umax=min(args.tail_u,1/(2*tt)) if tt>0 else args.tail_u
            extra.append(cyc(np.array(tt),np.array(1+v*(umax-1))))
    print({'dimension':n,'a':str(a),'b':str(b),'round':rnd,'raw_B':numeric['raw_B'],
           'guarded_B':numeric['guarded_B'],'extra_constraints':len(extra),'seconds':time.time()-start},flush=True)
    if not extra:break
    A=np.vstack([A,np.array(extra)]);rhs=np.r_[rhs,np.full(len(extra),-1e-3)]
