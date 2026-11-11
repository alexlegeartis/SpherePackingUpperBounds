if not __debug__:
    raise RuntimeError("Certificate verification requires assertions; do not run Python with -O or -OO.")

import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json,hashlib,time,sys
import numpy as np
from scipy.optimize import minimize,differential_evolution
ROOT=Path(__file__).resolve().parent; SRC=ROOT.parent/'density_kissing16_reweight_frozen_20261111'
assert hashlib.sha256((SRC/'candidate.json').read_bytes()).hexdigest()=='b3225ad4cce819395ce3bd47613f15697af4646fb722ded1916020b5e52ba573'
assert hashlib.sha256((SRC/'factors.json').read_bytes()).hexdigest()=='6860eca2c82fdb1b1ebb5d3a4108a216d46c900d6fac7bdeb0d880272e027954'
D=json.loads((SRC/'factors.json').read_text()); W=json.loads((SRC/'candidate.json').read_text())['record']; C=5003/10000
A=np.array(W['global_coefficients']); R=np.array(W['local_angular_coefficients']); T=np.array(W['theta'])
F=[]
for block in D['blocks']:
    L=np.zeros((block['size'],block['size']))
    for i,j,s in block['L_entries']:L[i,j]=float(s)
    F.append(L.T@L+float(block['lambda_decimal'])*np.eye(block['size']))
B=F[16]; F=F[17:]
def gegen(x,n,d):
    x=np.asarray(x); p=np.empty((*x.shape,d+1));p[...,0]=1
    if d:p[...,1]=x
    for k in range(1,d):p[...,k+1]=((2*k+n-2)*x*p[...,k]-k*p[...,k-1])/(k+n-2)
    return p

def kernel(u,v,t):
    u,v,t=np.broadcast_arrays(u,v,t);g=t-u*v;z=(1-u*u)*(1-v*v)
    q0=np.ones_like(u); q1=g;val=np.full_like(u,B[1,1])
    for k in range(17):
        if k==0:q=q0
        elif k==1:q=q1
        else:q0,q1=q1,((2*k+11)*g*q1-(k-1)*z*q0)/(k+12);q=q1
        if T[k]:
            p=gegen(u,16+2*k,16-k);r=gegen(v,16+2*k,16-k)
            val+=T[k]*q*np.einsum('...i,ij,...j->...',p,F[k],r,optimize=True)
    return val

def cycle(x):
    x=np.asarray(x);u,v,t=x.T
    return kernel(u,v,t)+kernel(u,t,v)+kernel(v,t,u)+(gegen(u,16,24)[...,1:]+gegen(v,16,24)[...,1:]+gegen(t,16,24)[...,1:])@R

def det(x):
    x=np.asarray(x);u,v,t=x.T;return 1+2*u*v*t-u*u-v*v-t*t

def feasible_x(rng,n,kind):
    u=rng.uniform(-1,C,n); v=u+(C-u)*rng.random(n)
    if kind=='cface':v[:]=C
    if kind=='equaluv':v=u.copy()
    rad=np.sqrt(np.maximum(0,(1-u*u)*(1-v*v)))
    lo=np.maximum(v,u*v-rad);hi=np.minimum(C,u*v+rad);valid=lo<=hi
    u,v,lo,hi=[x[valid] for x in (u,v,lo,hi)]
    if kind=='detlow':t=lo
    elif kind=='dethigh':t=hi
    elif kind=='cface':t=hi
    elif kind=='equaluv':t=lo+(hi-lo)*rng.random(len(u))
    elif kind=='equalvt':t=lo
    elif kind=='nearboundary':
        r=rng.random(len(u));r=np.where(r<.5,r**8,1-(1-r)**8);t=lo+(hi-lo)*r
    else:t=lo+(hi-lo)*rng.random(len(u))
    return np.c_[u,v,t]

def challenge():
    start=time.monotonic();rng=np.random.default_rng(74169282);best=[];count=0;records=[]
    for kind in ['random','nearboundary','detlow','dethigh','cface','equaluv','equalvt']:
        for batch in range(4):
            x=feasible_x(rng,65000,kind);vals=cycle(x);count+=len(x);ids=np.argsort(vals)[-20:];best.extend([(float(vals[i]),x[i]) for i in ids]);best=sorted(best,key=lambda z:z[0])[-150:]
        rec={'kind':kind,'evaluated_total':count,'maximum':best[-1][0],'point':best[-1][1].tolist(),'gram_det':float(det(best[-1][1])),'seconds':time.monotonic()-start};records.append(rec);print(rec,flush=True)
    # distinct starts: random high values and geometry edge structure
    starts=[]
    for val,x in sorted(best,key=lambda z:-z[0]):
        if all(np.linalg.norm(x-y)>.01 for y in starts):starts.append(x)
    starts+=list(feasible_x(rng,200,'random')[:40])
    for i,x in enumerate(starts):
        opt=minimize(lambda z:-float(cycle(z)),x,method='SLSQP',bounds=[(-1,C)]*3,constraints=[{'type':'ineq','fun':det}],options={'maxiter':250,'ftol':1e-13})
        if det(opt.x)>=-1e-11 and np.max(opt.x)<=C+1e-14:
            best.append((-float(opt.fun),opt.x));best=sorted(best,key=lambda z:z[0])[-150:]
        if (i+1)%10==0:print({'local_completed':i+1,'max':best[-1][0],'point':best[-1][1].tolist(),'det':float(det(best[-1][1]))},flush=True)
    # Differential evolution on full box, penalty for infeasible triples.
    def penalized(x):
        v=float(cycle(x));g=det(x)
        return -v+1e5*max(0,-g)
    for seed in range(3):
        de=differential_evolution(penalized,[(-1,C)]*3,popsize=12,maxiter=220,tol=1e-10,seed=seed+96491,polish=False)
        opt=minimize(lambda z:-float(cycle(z)),de.x,method='SLSQP',bounds=[(-1,C)]*3,constraints=[{'type':'ineq','fun':det}],options={'maxiter':250,'ftol':1e-13})
        if det(opt.x)>=-1e-11:best.append((-float(opt.fun),opt.x));best=sorted(best,key=lambda z:z[0])[-150:]
        print({'de_seed':seed,'localmax':-float(opt.fun),'point':opt.x.tolist(),'det':float(det(opt.x))},flush=True)
    output={'scope':'Independent numerical challenge only. No continuum proof. Optimization points with tiny negative determinant are not proof-acceptance points.','random_seed':74169282,'count':count,'sampling':records,'best':[{'value':v,'point':x.tolist(),'gram_det':float(det(x))} for v,x in best[-30:]],'seconds':time.monotonic()-start}
    (ROOT/'challenge.json').write_text(json.dumps(output,indent=2));print('FINISHED',output['seconds'],flush=True)
if __name__=='__main__':challenge()
