"""Sampled dimension-specific radial LP around the frozen source. NOT PROOF."""
from pathlib import Path
from fractions import Fraction as F
import mpmath as mp
import numpy as np,json,math,argparse
from scipy.optimize import linprog
from scipy.special import eval_genlaguerre,gammaln
ROOT=Path(__file__).resolve().parent;D=ROOT/'work';D.mkdir(exist_ok=True);S=ROOT/'inputs/source44';mp.mp.dps=160
ap=argparse.ArgumentParser();ap.add_argument('--dimension',type=int,default=43);ap.add_argument('--degree',type=int,default=31);ap.add_argument('--r2',type=float,default=5.5);args=ap.parse_args();n=args.dimension;d=args.degree;r2=args.r2;kappa=12.75/r2;R2=7;tag=str(r2).replace('.','p')
r=json.loads((S/'certificate.json').read_text());q=json.loads((S/'continuum_three_point/cyclic_angular40_exact.json').read_text());B=F(q['B_exact']);N48=F(r['fourier_polynomial'][0]);ps=[F(x)/N48 for x in r['physical_polynomial']];hs=[F(x)/N48 for x in r['fourier_polynomial']];pc=list(map(F,q['polynomial_coefficients']))
rf=lambda a,k:math.prod((a+j for j in range(k)),start=F(1))
HN=[sum((hs[k]*math.comb(k,j)*rf(F(44-n,2),k-j) for k in range(j,len(ps))),F(0)) for j in range(len(ps))]
b0=sum((pc[2*j]*rf(F(1,2),j)/rf(F(n,2),j) for j in range(len(pc)//2)),F(0));Bn=(float(F(json.loads((ROOT/'inputs/angular_exact.json').read_text())['B_exact'])) if True else float(F(json.loads((D/f'dim{n}/cyclic_angular40_exact.json').read_text())['B_exact'])) if (D/f'dim{n}/cyclic_angular40_exact.json').exists() else json.loads((D/f'angular{n}_a40_scaled_candidate.json').read_text())['B']*1.00001 if (D/f'angular{n}_a40_scaled_candidate.json').exists() else float(B/b0))
frac=lambda x:mp.mpf(x.numerator)/x.denominator
p=[frac(x) for x in ps];h=[frac(x) for x in HN]
norm=np.exp(gammaln(np.arange(d+1)+n/2)-gammaln(n/2)-gammaln(np.arange(d+1)+1));sign=(-1.)**np.arange(d+1)
def lag(x):return np.array([eval_genlaguerre(k,n/2-1,2*x)/norm[k] for k in range(d+1)])
def pol(cs,x):return float(mp.polyval(list(reversed(cs)),mp.mpf(float(x))))
cache=D/f'perturb_rows_{n}_{d}_r{tag}_B{Bn:.8g}.npz'
if cache.exists():
 A=np.load(cache)['A']
else:
 rows=[]
 for X in np.linspace(r2*np.pi,R2*np.pi,1001):
  hv=np.exp(-kappa*(X-r2*np.pi));g=np.exp(X)*(2*hv-hv*hv)/Bn
  rows.append(np.r_[pol(p,X),lag(X),-g])
 for X in np.r_[np.linspace(R2*np.pi,200,2001),np.geomspace(200,2000,401)]:rows.append(np.r_[pol(p,X),lag(X),0])
 for X in np.r_[np.linspace(0,200,2401),np.geomspace(200,2000,401)]:rows.append(np.r_[-pol(h,X),-lag(X)*sign,0])
 A=np.array(rows);scale=np.maximum(np.max(np.abs(A),axis=1),1e-100);A/=scale[:,None]
 np.savez_compressed(cache,A=A)
extra=D/f'radial_scaled{n}_degree{d}_r{tag}_extra.json';reject=D/f'radial_scaled{n}_degree{d}_r{tag}_rejection.json'
points=json.loads(extra.read_text()) if extra.exists() else {'physical':[],'fourier':[]}
if reject.exists():
 bad=json.loads(reject.read_text())['global_sign_roots']
 for name,roots in bad.items():
  roots=sorted(z[0] for z in roots)
  points[name]+=roots+[(x+y)/2 for x,y in zip(roots,roots[1:])]
 extra.write_text(json.dumps(points,indent=2))
rows=[]
for X in points['physical']:rows.append(np.r_[pol(p,X),lag(X),0])
for X in points['fourier']:rows.append(np.r_[-pol(h,X),-lag(X)*sign,0])
if rows:
 Ae=np.array(rows);Ae/=np.maximum(np.max(np.abs(Ae),axis=1,keepdims=True),1e-100);A=np.vstack([A,Ae])
c=np.r_[float(ps[0]),np.ones(d+1),1];eq=np.r_[float(HN[0]),sign,0]
res=linprog(c,A_ub=A,b_ub=np.full(len(A),-1e-8),A_eq=eq[None,:],b_eq=[1],bounds=[(0,None)]+[(None,None)]*(d+1)+[(0,None)],method='highs',options={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9})
if not res.success:print(res.message);raise SystemExit(1)
x=res.x;dens=res.fun*math.pi**(n/2)*(r2/4)**(n/2)/math.gamma(n/2+1)
out={'r2':r2,'R2':R2,'kappa':kappa,'dimension':n,'degree':d,'certified':False,'density_candidate':dens,'objective':res.fun,'x':x.tolist(),'B_n':Bn,'source_t':float(x[0]),'origin_cost':float(x[-1]),'max_violation':float(np.max(A@x))};print(out,flush=True)
(D/f'radial_scaled{n}_degree{d}_r{tag}.json').write_text(json.dumps(out,indent=2))
