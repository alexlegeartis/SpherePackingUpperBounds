"""Exact-polynomial far replay with independently certified stepwise tail budgets.
Only the far condition is checked here. No density claim follows by itself.
"""
if not __debug__:raise RuntimeError('Do not disable assertions.')
from pathlib import Path
from fractions import Fraction as F
import argparse,json,math,time,sys,importlib.util
import numpy as np
from flint import fmpq as Q,fmpq_poly as P,fmpq_mat as M,arb,ctx
ctx.prec=256
ap=argparse.ArgumentParser();ap.add_argument('witness',type=Path);ap.add_argument('budget',type=Path);ap.add_argument('--limit',type=int,default=100000);args=ap.parse_args()
data=json.loads(args.witness.read_text());src=json.loads((args.witness.parent/data['source']).read_text());n=src['dimension'];p=src['p'];K=src['K'];R=F(src['R']);assert R==F(3,2) and src['J']==1
budget=json.loads(args.budget.read_text());assert budget['schema']=='tail_funded_lens_step_allowances_v1'
records=[z for z in budget['dimensions'] if z['dimension']==n];assert len(records)==1;budget=records[0]
assert n==5 and p==4 and K==12 and src['kernel_family']=='tail_funded_lens_cyclic'
assert budget['dimension']==n and budget['geometry_mode']=='spherical_transverse'
expected_payload=dict(dimension=n,basis='L_k^(n/2-1)(2X)/L_k^(n/2-1)(0)',variable='X=pi*r^2',physical_factor='exp(-X)',additional_physical_terms='none',degree=31,exact_rational_coefficients=[str(F(str(c))) for c in src['x']])
assert budget['physical_function']==expected_payload
assert len(src['x'])==32 and F(budget['width'])==F(1,200) and len(budget['bins'])==132
assert [(F(row['lower_base']),row['cap'])for row in budget['exact_cap_geometry']]==[(F(3,2),10),(F(31,20),9),(F(8,5),8),(F(41,25),5),(F(17,10),5),(F(7,4),5),(F(9,5),3),(F(46,25),2),(F(48,25),1)]
assert list(map(F,budget['physical_function']['exact_rational_coefficients']))==list(map(lambda z:F(str(z)),src['x']))
assert budget['physical_function']['additional_physical_terms']=='none'
bins=[(F(z['lo']),F(z['hi']),int(z['cap']),F(z['eta']))for z in budget['bins']]
bins.sort();assert all(a<b and b-a==F(1,200) and eta>=0 for a,b,B,eta in bins)
assert bins[0][0]==F(3,2) and bins[-1][1]==F(54,25)
assert all(bins[i][1]==bins[i+1][0]for i in range(len(bins)-1))
rmax=F(27,25)
# Recompute every embedded geometry/code certificate and 1D budget gate.
# The research helper is explicitly adjacent to the supplied budget artifact.
verifier=Path(__file__).resolve().parent/'construct_step_allowances.py'
sys.path.insert(0,str(verifier.parent))
spec=importlib.util.spec_from_file_location('lens_budget_replay',verifier);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
assert module.PhysicalFunction(payload=expected_payload).physical_sha256==budget['physical_function_sha256']
assert module.validate_cap_geometry(budget)==budget['exact_cap_geometry']
assert module.replay(args.budget)['verified']
assert F(budget['short_radius'])==rmax
for a,b,B,eta in bins:assert a>=R and b<=2*rmax and B>=1
Qs=[M([[Q(z)for z in row]for row in A])for A in data['matrices']]
assert len(Qs)==K+1 and all(A==A.transpose()for A in Qs)
X=P([0,1]);Ts=[P([1]),2*X-1]
for i in range(1,p):Ts.append(2*Ts[1]*Ts[-1]-Ts[-2])
top=1-F(1,2)/(R*R);T=-1+Q(str(top+1))*X;G=[P([1]),T]
for k in range(1,K):G.append(((2*k+n-2)*T*G[-1]-k*G[-2])/Q(k+n-2))
def bern(q,d):return [sum((q[j]*Q(math.comb(i,j),math.comb(d,j))for j in range(min(i,q.degree())+1)),Q(0))for i in range(d+1)]
Ba=[bern(G[k],K)for k in range(K+1)]
def init(lo1,hi1,lo2,hi2):
 A=M([bern(q(Q(str(lo1))+Q(str(hi1-lo1))*X),p)for q in Ts]);B=M([bern(q(Q(str(lo2))+Q(str(hi2-lo2))*X),p)for q in Ts])
 blocks=[]
 for k,Qk in enumerate(Qs):
  o=int(k==0);qr=M([[Qk[i+o,j+o]for j in range(p+1)]for i in range(p+1)]);blocks.append(A.transpose()*qr*B)
 lo=np.empty((p+1,p+1,K+1));hi=lo.copy()
 for i,j,l in np.ndindex(lo.shape):
  v=arb(sum((blocks[k][i,j]*Ba[k][l]for k in range(K+1)),Q(0)))
  lo[i,j,l]=np.nextafter(float(v.lower()),-np.inf);hi[i,j,l]=np.nextafter(float(v.upper()),np.inf)
 return lo,hi

def split(l,h,axis):
 a=np.moveaxis(l,axis,0).copy();b=np.moveaxis(h,axis,0).copy();d=a.shape[0]-1;ll=np.empty_like(a);lh=np.empty_like(b);rl=np.empty_like(a);rh=np.empty_like(b);ll[0]=a[0];lh[0]=b[0];rl[d]=a[-1];rh[d]=b[-1]
 for k in range(1,d+1):
  a=np.nextafter(np.nextafter(a[:-1]+a[1:],-np.inf)*.5,-np.inf);b=np.nextafter(np.nextafter(b[:-1]+b[1:],np.inf)*.5,np.inf);ll[k]=a[0];lh[k]=b[0];rl[d-k]=a[-1];rh[d-k]=b[-1]
 return tuple(np.moveaxis(z,0,axis)for z in [ll,lh,rl,rh])

def distance_range(box):
 rs=[tuple(1+(R-1)*q for q in box[i])for i in [0,1]];tl,th=[-1+(top+1)*q for q in box[2]]
 def w(r,s,t):return r*r+s*s-2*r*s*t
 whi=max(w(r,s,tl)for r in rs[0]for s in rs[1]);candidates=[w(r,s,th)for r in rs[0]for s in rs[1]]
 def clip(x,a,b):return max(a,min(b,x))
 for r in rs[0]:candidates.append(w(r,clip(th*r,*rs[1]),th))
 for s in rs[1]:candidates.append(w(clip(th*s,*rs[0]),s,th))
 return min(candidates),whi,rs

def allowance(wlo,whi,rs):
 if rs[0][1]>rmax or rs[1][1]>rmax:return F(0)
 if wlo<=bins[0][0]**2 or whi>=bins[-1][1]**2:return F(0)
 return min(eta for a,b,B,eta in bins if a*a<=whi and b*b>=wlo)

cut=(rmax-1)/(R-1);intervals=[(F(0),cut),(cut,F(1))];stack=[]
for i,(a,b)in enumerate(intervals):
 for c,d in intervals[i:]:
  lo,hi=init(a,b,c,d);stack.append((lo,hi,((a,b),(c,d),(F(0),F(1))),0))
checked=passed=relaxed=excluded=0;maxdepth=0;fail=None;start=time.time()
while stack:
 lo,hi,box,depth=stack.pop();checked+=1;maxdepth=max(maxdepth,depth)
 if np.max(hi)<0:passed+=1;continue
 wlo,whi,rs=distance_range(box)
 if whi<R*R:excluded+=1;continue
 eta=allowance(wlo,whi,rs)
 if eta and np.max(hi)<np.nextafter(float(arb(Q(str(eta))).lower()),-np.inf):passed+=1;relaxed+=1;continue
 if checked>=args.limit or depth>=65:
  fail={'box':box,'depth':depth,'upper':float(np.max(hi)),'eta':str(eta),'wlo':str(wlo),'whi':str(whi),'remaining':len(stack)};break
 mid=(lo+hi)*.5;scores=[float(np.max(np.abs(np.diff(mid,axis=i))))for i in range(3)];axis=int(np.argmax(scores));ll,lh,rl,rh=split(lo,hi,axis);a,b=box[axis];m=(a+b)/2;bx1=list(box);bx2=list(box);bx1[axis]=(a,m);bx2[axis]=(m,b);stack.extend([(rl,rh,tuple(bx2),depth+1),(ll,lh,tuple(bx1),depth+1)])
 if checked%1000==0:print('PROGRESS',checked,'relaxed',relaxed,'stack',len(stack),flush=True)
out={'verified':fail is None and not stack,'checked':checked,'passed':passed,'relaxed_boxes':relaxed,'excluded':excluded,'max_depth':maxdepth,'failure':fail,'seconds':time.time()-start,'witness':str(args.witness),'budget':str(args.budget),'scope':'Far domain only; step budgets must also be replayed and all other density hypotheses verified.'}
print(json.dumps(out,default=str),flush=True);args.witness.with_name('witness_wedge_far.json').write_text(json.dumps(out,indent=2,default=str)+'\n')
