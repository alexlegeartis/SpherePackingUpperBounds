"""Numerical radial continuation of the verified44D Gaussian cyclic kernel.
All finite samples are discovery only. Exact global verification is required.
"""
import argparse,json,math,sys
from pathlib import Path
from fractions import Fraction as F
import numpy as np,mpmath as mp
from scipy.optimize import linprog
from scipy.special import eval_genlaguerre,gammaln
sys.set_int_max_str_digits(0);mp.mp.dps=140
D=Path(__file__).resolve().parent;ROOT=D.parent
ap=argparse.ArgumentParser();ap.add_argument('--dimension',type=int,default=43);ap.add_argument('--degree',type=int,default=31);ap.add_argument('--r2',type=str,default='23/4');ap.add_argument('--angular-candidate');args=ap.parse_args();n=args.dimension;d=args.degree;r2=F(args.r2);kappa=F(51,4)/r2
r=json.loads((D/'source_inputs/source44.json').read_text());q=json.loads((D/'source_inputs/angular44.json').read_text())
rf=lambda a,k:math.prod((a+j for j in range(k)),start=F(1))
norm0=F(r['fourier_polynomial'][0]);ps=[F(x)/norm0 for x in r['physical_polynomial']];hs=[F(x)/norm0 for x in r['fourier_polynomial']];pc=list(map(F,q['polynomial_coefficients']))
HN=[sum((hs[k]*math.comb(k,j)*rf(F(44-n,2),k-j) for k in range(j,len(ps))),F(0)) for j in range(len(ps))]
b0=sum((pc[2*j]*rf(F(1,2),j)/rf(F(n,2),j) for j in range(len(pc)//2)),F(0));Bn=F(q['B_exact'])/b0
if args.angular_candidate: Bn=F(str(json.loads(Path(args.angular_candidate).read_text())['B']))*F(1002,1000)
frac=lambda v:mp.mpf(v.numerator)/v.denominator
p=list(map(frac,ps));h=list(map(frac,HN));norm=np.exp(gammaln(np.arange(d+1)+n/2)-gammaln(n/2)-gammaln(np.arange(d+1)+1));sign=(-1.)**np.arange(d+1)
def lag(x):return np.array([eval_genlaguerre(k,n/2-1,2*x)/norm[k] for k in range(d+1)])
def pol(cs,x):return float(mp.polyval(list(reversed(cs)),mp.mpf(float(x))))
rows=[]
for X in np.linspace(float(r2)*np.pi,7*np.pi,1601):
 hv=np.exp(-float(kappa)*(X-float(r2)*np.pi));g=np.exp(X)*(2*hv-hv*hv)/float(Bn)
 rows.append(np.r_[pol(p,X),lag(X),-g])
for X in np.r_[np.linspace(7*np.pi,250,4001),np.geomspace(250,2500,401)]:rows.append(np.r_[pol(p,X),lag(X),0])
for X in np.r_[np.linspace(0,250,4001),np.geomspace(250,2500,401)]:rows.append(np.r_[-pol(h,X),-lag(X)*sign,0])
tag=f'radial{n}_d{d}_r{str(r2).replace("/","_")}{"_newangular" if args.angular_candidate else ""}'
extra=D/(tag+'_extra.json');reject=D/(tag+'_rejection.json')
points=json.loads(extra.read_text()) if extra.exists() else {'physical':[],'fourier':[]}
if reject.exists():
 for name,roots in json.loads(reject.read_text())['global_sign_roots'].items():
  zs=sorted(v[0] for v in roots);points[name]+=zs+[(x+y)/2 for x,y in zip(zs,zs[1:])]
 extra.write_text(json.dumps(points,indent=2))
for X in points['physical']:rows.append(np.r_[pol(p,X),lag(X),0])
for X in points['fourier']:rows.append(np.r_[-pol(h,X),-lag(X)*sign,0])
A=np.array(rows);A/=np.maximum(np.max(abs(A),axis=1,keepdims=True),1e-100)
c=np.r_[float(ps[0]),np.ones(d+1),1];eq=np.r_[float(HN[0]),sign,0]
res=linprog(c,A_ub=A,b_ub=np.full(len(A),-1e-8),A_eq=eq[None,:],b_eq=[1],bounds=[(0,None)]+[(None,None)]*(d+1)+[(0,None)],method='highs',options={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9})
if not res.success:print(res.message);raise SystemExit(1)
dens=res.fun*math.pi**(n/2)*(float(r2)/4)**(n/2)/math.gamma(n/2+1)
out=dict(dimension=n,degree=d,r2=str(r2),R2='7',kappa=str(kappa),angular_B_exact=str(Bn),source='dim44_gaussian_ready_20261111',certified=False,density_candidate=dens,x=res.x.tolist(),objective=res.fun,max_violation=float(np.max(A@res.x)))
(D/f'radial{n}_d{d}_r{str(r2).replace("/","_")}{"_newangular" if args.angular_candidate else ""}.json').write_text(json.dumps(out,indent=2));print(out,flush=True)
