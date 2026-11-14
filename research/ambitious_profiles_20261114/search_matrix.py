"""Two-radial-feature PSD kernel SDP. Every output is discovery, never proof.

Unlike the scalar angular family, the finite-scale angular kernel may be
positive on part of the old negative cap. Only the slow asymptote and its
slow/fast endpoint combination retain cap constraints.
"""
import argparse,json,time,math,sys
from pathlib import Path
from fractions import Fraction as F
import numpy as np
import cvxpy as cp
from scipy.special import eval_genlaguerre,gammaln
from flint import arb,arb_poly,ctx
from common import HERE,base_angular,source

ap=argparse.ArgumentParser()
ap.add_argument('dimension',type=int);ap.add_argument('--degree',type=int,default=63)
ap.add_argument('--s',default='6');ap.add_argument('--slow',default='40');ap.add_argument('--fast',default='60')
ap.add_argument('--cut',type=float,default=40);ap.add_argument('--round',type=int,default=0)
ap.add_argument('--fixed-leading',action='store_true')
ap.add_argument('--fixed-radial')
args=ap.parse_args();n,d=args.dimension,args.degree;s=F(args.s);slow,fast=F(args.slow),F(args.fast)
assert 0<slow<fast
D=HERE/(f'matrix_n{n}_d{d}_s{str(s).replace("/","_")}_b{slow}_{fast}'+('_fixed' if args.fixed_leading else ''))
if args.fixed_radial:D=Path(str(D)+'_fr')
D.mkdir(exist_ok=True);H=D/f'round{args.round:02}';H.mkdir(exist_ok=True)
(H/'search_source.py').write_text(Path(__file__).read_text())
ctx.prec=512;start=time.time()
ref=json.loads(base_angular(n).read_text());bell=np.array([float(F(x)) for x in ref['gegenbauer_coefficients']])
B=float(F(ref['B_exact']));weights=bell/B;ellmax=len(bell)-1
active=[k for k,z in enumerate(bell) if z]
p,h,ls=source(n,d);AP=arb_poly([arb(x) for x in p.coeffs()]);AH=arb_poly([arb(x) for x in h.coeffs()])
norm=np.exp(gammaln(np.arange(d+1)+n/2)-gammaln(n/2)-gammaln(np.arange(d+1)+1));sign=(-1.)**np.arange(d+1)
def lag(x):return np.array([eval_genlaguerre(k,n/2-1,2*np.asarray(x))/norm[k] for k in range(d+1)]).T
def pol(poly,x):return np.array([float(poly(arb(float(z))).mid()) for z in x])
def gegen(t):
    t=np.asarray(t);g=[np.ones_like(t),t]
    for k in range(1,ellmax):g.append(((2*k+n-2)*t*g[-1]-k*g[-2])/(k+n-2))
    return np.stack(g,axis=-1)
def tri_rows(t,u,T):
    t,u,T=np.broadcast_arrays(np.asarray(t),np.asarray(u),np.asarray(T))
    w2=1+u*u-2*t*u;w=np.sqrt(w2)
    g0,g1,g2=gegen(t),gegen((1-u*t)/w),gegen((u-t)/w)
    w1=np.exp(-float(slow)*T*(w2-u*u));w2g=np.exp(-float(slow)*T*(w2-1))
    r0=np.exp(-float(fast-slow)*(T-1));r1=np.exp(-float(fast-slow)*(T*u*u-1));r2=np.exp(-float(fast-slow)*(T*w2-1))
    aa=g0+w1[...,None]*g1+w2g[...,None]*g2
    ab=(r0+r1)[...,None]*g0+(w1*(r0+r2))[...,None]*g1+(w2g*(r1+r2))[...,None]*g2
    bb=(r0*r1)[...,None]*g0+(w1*r0*r2)[...,None]*g1+(w2g*r1*r2)[...,None]*g2
    return np.stack([aa,ab,bb],axis=-1)[...,active,:].reshape(-1,3*len(active))*np.repeat(weights[active],3)
def leading_rows(t,u):
    w2=1+u*u-2*t*u;w=np.sqrt(w2)
    vals=gegen(t)+np.exp(-float(slow)*(w2-u*u))[:,None]*gegen((1-u*t)/w)+np.exp(-float(slow)*(w2-1))[:,None]*gegen((u-t)/w)
    z=np.zeros((len(t),len(active),3));z[:,:,0]=vals[:,active]
    return z.reshape(len(t),-1)*np.repeat(weights[active],3)

T0,U0=np.meshgrid(np.r_[np.linspace(-1,0,31),np.linspace(0,.5,51)],np.r_[np.linspace(1,1.3,21),np.linspace(1.3,2,13),np.linspace(2,5,9)])
ts,us=T0.ravel(),U0.ravel();mask=2*ts*us<=1;ts,us=ts[mask],us[mask]
scales=np.array([1,1.025,1.075,1.15,1.3,1.6,2,3])
tri=np.vstack([tri_rows(ts,us,np.full_like(ts,T)) for T in scales])
lead=leading_rows(ts,us)
cap=gegen(np.linspace(-1,.5001,2001))[:,active]*weights[active]
gridcap=np.zeros((len(cap),len(active),3));gridcap[:,:,0]=cap
cap0=gridcap.reshape(len(cap),-1).copy()
gridcap[:,:,1]=cap
capsum=gridcap.reshape(len(cap),-1).copy()
extra=D/'cuts.json'
cuts=json.loads(extra.read_text()) if extra.exists() else {}
if cuts.get('triangles'):
    zz=np.array(cuts['triangles']);tri=np.vstack([tri,tri_rows(zz[:,0],zz[:,1],zz[:,2])])
if cuts.get('leading'):
    zz=np.array(cuts['leading']);lead=np.vstack([lead,leading_rows(zz[:,0],zz[:,1])])
if cuts.get('cap'):
    zz=gegen(np.array(cuts['cap']))[:,active]*weights[active]
    g=np.zeros((len(zz),len(active),3));g[:,:,0]=zz
    cap0=np.vstack([cap0,g.reshape(len(zz),-1)]);g[:,:,1]=zz
    capsum=np.vstack([capsum,g.reshape(len(zz),-1)])
for matrix in [tri,lead,cap0,capsum]:
    matrix/=np.maximum(np.max(abs(matrix),axis=1,keepdims=True),1e-200)

X=cp.Variable((3,3),symmetric=True,name='origin_block')
A={0:X[1:,1:]}
for k in active[1:]:A[k]=cp.Variable((2,2),symmetric=True,name=f'A{k}')
kernel=cp.hstack([cp.hstack([A[k][0,0],A[k][0,1],A[k][1,1]]) for k in active])
if args.fixed_radial:
    fixed=json.loads(Path(args.fixed_radial).read_text())
    fixed_polynomial='physical_polynomial' in fixed
    if fixed_polynomial:
        fixed_H0=F(fixed['fourier_polynomial'][0]);assert fixed_H0>0
        fixed_coeffs=[F(z)/fixed_H0 for z in fixed['physical_polynomial']]
        fixed_AP=arb_poly([arb(z.numerator)/z.denominator for z in fixed_coeffs])
        fixed_f0=float(fixed_AP(arb(0)).mid());fixed_values=[0.]*(d+2)
    else:fixed_values=fixed.get('radial',fixed.get('x',[])[:-1])
    assert fixed['dimension']==n and len(fixed_values)<=d+2
    fixed_values=list(fixed_values)+[0.]*(d+2-len(fixed_values))
    radial=cp.Constant(np.array(fixed_values))
else:radial=cp.Variable(d+2,name='radial')
cons=[X>>0,X[0,1]<=-1e-5,A[0][0,0]>=1e-5]
if not args.fixed_radial:cons.append(radial[0]>=0)
cons.extend(A[k]>>0 for k in active[1:])
if args.fixed_leading:
    cons.extend(A[k][0,0]==A[0][0,0] for k in active[1:])
cons.extend([tri@kernel<=-1e-7,lead@kernel<=-1e-7,cap0@kernel<=-1e-7,capsum@kernel<=-1e-7])
limit=max(300.,4*d+n+20.)
fouriergrid=np.r_[np.linspace(0,limit,4001),np.geomspace(limit,2500,301),cuts.get('fourier',[])]
physicalgrid=np.r_[np.linspace(args.cut,limit,3001),np.geomspace(limit,2500,301),cuts.get('physical',[])]
pairgrid=np.unique(np.r_[np.linspace(float(s)*np.pi,args.cut,1201),
             float(s)*np.pi+np.geomspace(1e-7,2.,151),cuts.get('pair',[])])
Frows=np.c_[-pol(AH,fouriergrid),-lag(fouriergrid)*sign]
Prows=np.c_[pol(AP,physicalgrid),lag(physicalgrid)]
Frows/=np.maximum(abs(Frows).max(axis=1,keepdims=True),1e-200)
Prows/=np.maximum(abs(Prows).max(axis=1,keepdims=True),1e-200)
if not args.fixed_radial:
    cons.extend([Frows@radial<=-1e-8,Prows@radial<=-1e-8,
                 np.r_[float(h[0]),sign]@radial==1])
v=pairgrid/(float(s)*np.pi)
h0=np.exp(-float(slow)*(v-1));h1=np.exp(-float(fast)*(v-1))
prow=np.c_[pol(AP,pairgrid),lag(pairgrid)]*(B*np.exp(-pairgrid))[:,None]
diag=np.tile(np.c_[h0*h0,2*h0*h1,h1*h1],(1,len(active)))*np.repeat(weights[active],3)
scale=np.maximum(np.max(abs(prow),axis=1),np.maximum(np.max(abs(diag),axis=1),2*np.maximum(h0,h1)))
if args.fixed_radial and fixed_polynomial:
    fixed_pair=pol(fixed_AP,pairgrid)*B*np.exp(-pairgrid)
    scale=np.maximum(abs(fixed_pair),np.maximum(np.max(abs(diag),axis=1),2*np.maximum(h0,h1)))
    physical_expression=fixed_pair/scale
else:physical_expression=(prow/scale[:,None])@radial
cons.append(physical_expression+(diag/scale[:,None])@kernel+
            2*(h0/scale)*X[0,1]+2*(h1/scale)*X[0,2]<=-1e-9)
objective=(fixed_f0 if args.fixed_radial and fixed_polynomial else np.r_[float(p[0]),np.ones(d+1)]@radial)+X[0,0]
problem=cp.Problem(cp.Minimize(objective),cons)
solver_attempts=[]
for tolerance in [1e-10,1e-8,1e-7]:
    try:
        value=problem.solve(solver='CLARABEL',tol_feas=tolerance,tol_gap_abs=tolerance,tol_gap_rel=tolerance,
                            max_iter=200,max_threads=2,verbose=False)
        solver_attempts.append({'tolerance':tolerance,'status':problem.status})
        if X.value is not None:break
    except cp.error.SolverError as exc:
        solver_attempts.append({'tolerance':tolerance,'error':str(exc)})
(H/'solver_attempts.json').write_text(json.dumps(solver_attempts,indent=2)+'\n')
assert radial.value is not None and X.value is not None,problem.status
blocks=[np.zeros((2,2)).tolist() for _ in range(ellmax+1)]
for k in active:blocks[k]=A[k].value.tolist()
out={'dimension':n,'degree':d,'s_exact':str(s),'slow_exact':str(slow),'fast_exact':str(fast),
     'fixed_slow_angular':args.fixed_leading,
     'fixed_radial_input':args.fixed_radial,
     'cut_X_exact':str(F(str(args.cut))),'reference_angular':str(base_angular(n).relative_to(HERE.parents[1])),
     'radial':radial.value.tolist(),'origin_block':X.value.tolist(),'blocks':blocks,
     'objective':float(value),'density_candidate':float(value)*np.pi**(n/2)*(float(s)/4)**(n/2)/math.gamma(n/2+1),
     'certified':False,'solver_status':problem.status,'seconds':time.time()-start,
     'solver_attempts':solver_attempts,
     'minimum_eigenvalue':min(np.linalg.eigvalsh(X.value).min(),*(np.linalg.eigvalsh(A[k].value).min() for k in active[1:])),
     'max_triangle_residual':float(max(tri@kernel.value)),
     'max_cap_residual':float(max(max(cap0@kernel.value),max(capsum@kernel.value))),
     'origin_cost':float(X.value[0,0]),'slow_cross':float(X.value[0,1]),'fast_cross':float(X.value[0,2]),
     'triangle_rows':len(tri),'leading_rows':len(lead),'pair_rows':len(pairgrid)}
(H/'candidate.json').write_text(json.dumps(out,indent=2)+'\n')
(D/'candidate.json').write_text(json.dumps(out,indent=2)+'\n')
print({k:v for k,v in out.items() if k not in ['radial','origin_block','blocks']},flush=True)
