"""Finite LP discovery, exact-root cuts and numerical pair separation. NOT PROOF."""
import argparse, json, math, os, sys, time
from pathlib import Path
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog, minimize_scalar
from scipy.special import eval_genlaguerre, gammaln
from flint import arb_poly, arb, ctx
from model import inputs

ap=argparse.ArgumentParser()
ap.add_argument('dimension',type=int,choices=[45,46,47])
ap.add_argument('--degree',type=int,default=63)
ap.add_argument('--cut',type=float,default=40)
ap.add_argument('--margin',type=float,default=1e-8)
args=ap.parse_args()
n,d,cut=args.dimension,args.degree,args.cut
ctx.prec=512
D=Path(__file__).resolve().parent/f'n{n}_d{d}_c{cut:g}'
D.mkdir(exist_ok=True)
start=time.time()
p,h,ls,B=inputs(n,d)
apoly=arb_poly([arb(c) for c in p.coeffs()])
ahpoly=arb_poly([arb(c) for c in h.coeffs()])
norm=np.exp(gammaln(np.arange(d+1)+n/2)-gammaln(n/2)-gammaln(np.arange(d+1)+1))
sign=(-1.)**np.arange(d+1)
def lag(x):
    return np.array([eval_genlaguerre(k,n/2-1,2*x)/norm[k] for k in range(d+1)])
def pol(poly,x):
    v=poly(arb(float(x)))
    assert v.is_finite() and v.rel_accuracy_bits()>80, (x,str(v))
    return float(v.mid())
def row(kind,x):
    if kind=='fourier':
        return np.r_[-pol(ahpoly,x),-lag(x)*sign,0]
    if kind=='physical':
        return np.r_[pol(apoly,x),lag(x),0]
    hv=np.exp(-17/8*(x-6*np.pi))
    allowance=np.exp(x)*(2*hv-hv*hv)/float(B)
    return np.r_[pol(apoly,x),lag(x),-allowance]

cache=D/'rows.npz'
if cache.exists():
    data=np.load(cache); A=data['A']
else:
    grid={
      'pair':np.unique(np.r_[np.linspace(6*np.pi,7*np.pi,1001),np.linspace(7*np.pi,cut,2001)]),
      'physical':np.r_[np.linspace(cut,250,3001),np.geomspace(250,2500,401)],
      'fourier':np.r_[np.linspace(0,250,4001),np.geomspace(250,2500,401)]}
    A=np.array([row(kind,float(x)) for kind,grid1 in grid.items() for x in grid1])
    A/=np.maximum(np.max(abs(A),axis=1,keepdims=True),1e-100)
    np.savez_compressed(cache,A=A)
extra=D/'cuts.json'
points=json.loads(extra.read_text()) if extra.exists() else {'physical':[],'fourier':[],'pair':[]}
if extra.exists():
    rows=np.array([row(kind,x) for kind,zs in points.items() for x in zs])
    if len(rows):
        rows/=np.maximum(np.max(abs(rows),axis=1,keepdims=True),1e-100)
        A=np.vstack([A,rows])
c=np.r_[float(p[0]),np.ones(d+1),1]
eq=np.r_[float(h[0]),sign,0]
options={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9}
attempts=[]
for method,presolve in [('highs',True),('highs-ipm',True),('highs-ds',False)]:
    res=linprog(c,A_ub=A,b_ub=np.full(len(A),-args.margin),A_eq=eq[None,:],b_eq=[1],
                bounds=[(0,None)]+[(None,None)]*(d+1)+[(0,None)],method=method,
                options={**options,'presolve':presolve})
    attempts.append({'method':method,'presolve':presolve,'status':res.message})
    if res.success:break
if not res.success:
    raise RuntimeError(res.message)
dens=res.fun*math.pi**(n/2)*1.5**(n/2)/math.gamma(n/2+1)
out={'dimension':n,'degree':d,'coupled_upper_X_exact':str(Fraction(str(cut))),'certified':False,
     'density_candidate':dens,'objective':res.fun,'x':res.x.tolist(),
     'B_exact':str(B),'margin':args.margin,'max_scaled_violation':float(np.max(A@res.x)),
     'seconds':time.time()-start,'python':sys.version,'solver_attempts':attempts}
(D/'candidate.json').write_text(json.dumps(out,indent=2)+'\n')
print({k:v for k,v in out.items() if k not in ['x','B_exact']},flush=True)
