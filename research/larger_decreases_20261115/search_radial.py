"""Large radial LP: scaled seed plus degree95/111+ corrections and free scale.
Every result is a numerical candidate, never a density bound.
"""
import argparse,json,time,math
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog
from scipy.special import eval_genlaguerre,gammaln
from flint import arb,arb_poly,ctx
from common import HERE,base_angular,source

ap=argparse.ArgumentParser()
ap.add_argument('dimension',type=int);ap.add_argument('--degree',type=int,default=95)
ap.add_argument('--s',default='6');ap.add_argument('--a',default='0');ap.add_argument('--b',default='40')
ap.add_argument('--cut',type=float,default=60);ap.add_argument('--angular')
ap.add_argument('--margin',type=float,default=1e-7);ap.add_argument('--tag',default='')
ap.add_argument('--q',default='1')
args=ap.parse_args();n,d=args.dimension,args.degree;s=F(args.s);a,b=F(args.a),F(args.b)
amplitude=F(args.q);assert amplitude>0
assert s>0 and a+b>0 and b>0 and args.cut>float(s)*np.pi
ang=Path(args.angular).resolve() if args.angular else base_angular(n)
angular=json.loads(ang.read_text());B=F(angular['B_exact'])
assert angular['dimension']==n
if 'profile_a_exact' in angular:
    assert a==F(angular['profile_a_exact']) and b==F(angular['profile_b_exact'])
else:assert a==0 and b>=40
tag=args.tag or f'n{n}_d{d}_s{str(s).replace("/","_")}_a{str(a).replace("/","_")}_b{str(b).replace("/","_")}_q{str(amplitude).replace("/","_")}_c{args.cut:g}'
D=HERE/tag;D.mkdir(exist_ok=True)
ctx.prec=768;start=time.time();p,h,ls=source(n,d)
apoly=arb_poly([arb(c) for c in p.coeffs()]);ahpoly=arb_poly([arb(c) for c in h.coeffs()])
norm=np.exp(gammaln(np.arange(d+1)+n/2)-gammaln(n/2)-gammaln(np.arange(d+1)+1));sign=(-1.)**np.arange(d+1)
def lag(x):
    return np.array([eval_genlaguerre(k,n/2-1,2*np.asarray(x))/norm[k] for k in range(d+1)]).T
def pol(poly,x):
    val=poly(arb(float(x)));assert val.is_finite() and val.rel_accuracy_bits()>80,(x,str(val))
    return float(val.mid())
def rows(kind,x):
    x=np.asarray(x);lagx=lag(x)
    if kind=='fourier':return np.c_[-np.array([pol(ahpoly,z) for z in x]),-lagx*sign,np.zeros(len(x))]
    if kind=='physical':return np.c_[np.array([pol(apoly,z) for z in x]),lagx,np.zeros(len(x))]
    v=x/(np.pi*float(s));hv=np.exp(-float(a)*np.log(v)-float(b)*(v-1))
    hv*=float(amplitude)
    allowance=np.exp(x)*(2*hv-hv*hv)/float(B)
    return np.c_[np.array([pol(apoly,z) for z in x]),lagx,-allowance]

limit=max(300.,4*d+n+20.)
grids={'pair':np.unique(np.r_[np.linspace(float(s)*np.pi,args.cut,2501),
                float(s)*np.pi+np.geomspace(1e-8,min(3.,args.cut-float(s)*np.pi),401)]),
       'physical':np.r_[np.linspace(args.cut,limit,4001),np.geomspace(limit,2500,401)],
       'fourier':np.r_[np.linspace(0,limit,5001),np.geomspace(limit,2500,401)]}
extra=D/'cuts.json'
if extra.exists():
    points=json.loads(extra.read_text())
    for kind,xx in points.items():grids[kind]=np.r_[grids[kind],xx]
A=np.vstack([rows(kind,x) for kind,x in grids.items()])
if d>=111:A[:,0]=0
A/=np.maximum(abs(A).max(axis=1,keepdims=True),1e-200)
objective=np.r_[float(p[0]),np.ones(d+1),1]
eq=np.r_[float(h[0]),sign,0]
attempts=[]
for method,presolve in [('highs',True),('highs-ipm',True),('highs-ds',False)]:
    bounds=[(0,None) if d<111 else (0,0)]+[(None,None)]*(d+1)+[(0,None)]
    if d>=111 and d%2:bounds[d+1]=(0,None)
    r=linprog(objective,A_ub=A,b_ub=np.full(len(A),-args.margin),A_eq=eq[None,:],b_eq=[1],
        bounds=bounds,method=method,
        options={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9,'presolve':presolve,'time_limit':120})
    attempts.append({'method':method,'presolve':presolve,'status':r.message})
    if r.success:break
if not r.success:raise RuntimeError(r.message)
density=r.fun*np.pi**(n/2)*(float(s)/4)**(n/2)/math.gamma(n/2+1)
out={'dimension':n,'degree':d,'s_exact':str(s),'a_exact':str(a),'b_exact':str(b),
     'q_exact':str(amplitude),
     'cut_X_exact':str(F(str(args.cut))),'B_exact':str(B),'angular_file':str(ang.relative_to(HERE.parent.parent)),
     'certified':False,'x':r.x.tolist(),'objective':float(r.fun),'density_candidate':float(density),
     'max_scaled_residual':float(max(A@r.x)),'margin':args.margin,'seconds':time.time()-start,'solver_attempts':attempts}
(D/'candidate.json').write_text(json.dumps(out,indent=2)+'\n')
print({k:v for k,v in out.items() if k not in ['x','B_exact']},flush=True)
