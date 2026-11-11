"""Sampled dimension-specific radial LP around the frozen source. NOT PROOF."""
from pathlib import Path
from fractions import Fraction as F
import mpmath as mp
import numpy as np,json,math,argparse
from scipy.optimize import linprog
from scipy.special import eval_genlaguerre,gammaln
D=Path(__file__).parent;S=D/'inputs';mp.mp.dps=160
ap=argparse.ArgumentParser();ap.add_argument('--dimension',type=int,default=46);ap.add_argument('--degree',type=int,default=15);ap.add_argument('--control-safe-margin',action='store_true');args=ap.parse_args();n=args.dimension;d=args.degree
r=json.loads((S/'adaptive_cyclic40_radial111.json').read_text());q=json.loads((D/f'cyclic_angular{n}_exact.json').read_text());B=F(q['B_exact']);N48=F(r['fourier_polynomial'][0]);ps=[F(x)/N48 for x in r['physical_polynomial']];hs=[F(x)/N48 for x in r['fourier_polynomial']];pc=list(map(F,q['polynomial_coefficients']))
rf=lambda a,k:math.prod((a+j for j in range(k)),start=F(1))
HN=[sum((hs[k]*math.comb(k,j)*rf(F(48-n,2),k-j) for k in range(j,len(ps))),F(0)) for j in range(len(ps))]
b0=sum((pc[2*j]*rf(F(1,2),j)/rf(F(n,2),j) for j in range(len(pc)//2)),F(0));Bn=float(B/b0)
frac=lambda x:mp.mpf(x.numerator)/x.denominator
p=[frac(x) for x in ps];h=[frac(x) for x in HN]
norm=np.exp(gammaln(np.arange(d+1)+n/2)-gammaln(n/2)-gammaln(np.arange(d+1)+1));sign=(-1.)**np.arange(d+1)
def lag(x):return np.array([eval_genlaguerre(k,n/2-1,2*x)/norm[k] for k in range(d+1)])
def pol(cs,x):return float(mp.polyval(list(reversed(cs)),mp.mpf(float(x))))
def row_scale(A):return np.maximum(np.maximum(np.max(abs(A[:,:min(d,31)+2]),axis=1),abs(A[:,-1])),1e-100)
cache=D/f'perturb_rows_{n}_{d}_ref31.npz'
if cache.exists():
 A=np.load(cache)['A']
else:
 rows=[]
 for X in np.linspace(6*np.pi,7*np.pi,1001):
  hv=np.exp(-17/8*(X-6*np.pi));g=np.exp(X)*(2*hv-hv*hv)/Bn
  rows.append(np.r_[pol(p,X),lag(X),-g])
 for X in np.r_[np.linspace(7*np.pi,200,2001),np.geomspace(200,2000,401)]:rows.append(np.r_[pol(p,X),lag(X),0])
 for X in np.r_[np.linspace(0,200,2401),np.geomspace(200,2000,401)]:rows.append(np.r_[-pol(h,X),-lag(X)*sign,0])
 A=np.array(rows);scale=row_scale(A);A/=scale[:,None]
 np.savez_compressed(cache,A=A)
extra=D/f'radial_perturb{n}_degree{d}_extra.json';reject=D/f'radial_perturb{n}_degree{d}_rejection.json'
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
 Ae=np.array(rows);Ae/=row_scale(Ae)[:,None];A=np.vstack([A,Ae])
c=np.r_[float(ps[0]),np.ones(d+1),1];eq=np.r_[float(HN[0]),sign,0]
baseline=json.loads((D/f'baseline{n}_degree31.json').read_text());old=np.array(baseline['x']);embedded=np.r_[old[0],old[1:33],np.zeros(d-31),old[-1]]
control={'dimension':n,'enriched_degree':d,'reference_row_scaling_degree':31,'padded_control_max_lhs':float(np.max(A@embedded)),'padded_control_normalization_residual':float(eq@embedded-1),'padded_control_objective':float(c@embedded),'saved_control_objective':baseline['objective'],'same_objective':bool(abs(c@embedded-baseline['objective'])<1e-11),'numerical_margin_pass':bool(np.max(A@embedded)<=-1e-7+2e-9)}
(D/f'control_embedding{n}.json').write_text(json.dumps(control,indent=2)+'\n');print('CONTROL_EMBEDDING',control,flush=True)
rhs=np.full(len(A),-1e-7)
if args.control_safe_margin:
 rhs=np.maximum(rhs,A@embedded+1e-12);assert np.max(rhs)<0
 control['control_safe_margin']=True;control['relaxed_margin_rows']=int(np.sum(rhs>-1e-7));control['least_negative_rhs']=float(np.max(rhs));control['control_max_constraint_residual']=float(np.max(A@embedded-rhs));(D/f'control_embedding{n}.json').write_text(json.dumps(control,indent=2)+'\n');print('CONTROL_SAFE_MARGIN',control,flush=True)
res=linprog(c,A_ub=A,b_ub=rhs,A_eq=eq[None,:],b_eq=[1],bounds=[(0,None)]+[(None,None)]*(d+1)+[(0,None)],method='highs',options={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9,'time_limit':90})
if not res.success:print(res.message);raise SystemExit(1)
x=res.x;dens=res.fun*math.pi**(n/2)*1.5**(n/2)/math.gamma(n/2+1)
out={'dimension':n,'degree':d,'certified':False,'density_candidate':dens,'objective':res.fun,'x':x.tolist(),'B_n':Bn,'source_t':float(x[0]),'origin_cost':float(x[-1]),'max_violation':float(np.max(A@x))};print(out,flush=True)
(D/f'radial_perturb{n}_degree{d}.json').write_text(json.dumps(out,indent=2))
