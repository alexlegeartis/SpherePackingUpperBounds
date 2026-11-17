"""Fresh numerical SDP reconstruction; output is NOT a certified bound."""
import json,sys,time,math
from pathlib import Path
from fractions import Fraction as F
import numpy as np
import cvxpy as cp
from cap_extrema import extrema
from scipy.special import eval_gegenbauer,eval_genlaguerre,gammaln
from scipy.optimize import differential_evolution,minimize_scalar
import itertools
from flint import arb,arb_poly,ctx,fmpq
sys.set_int_max_str_digits(0);ctx.prec=512
D=Path(__file__).parent;p=int(sys.argv[1])if len(sys.argv)>1 else 6;seed=51
v=json.load(open(D/'base/geometry/dim51/frozen.json'));n=51;q=float(F(v['q']));co=np.array([float(F(x))for x in v['co']]);weights=co/q;N=len(co);k=np.arange(N);den=eval_gegenbauer(k,(n-2)/2,1);ij=[(i,j)for i in range(3)for j in range(i,3)]
def gegen(t):return eval_gegenbauer(k,(n-2)/2,np.asarray(t)[...,None])/den

def bilrows(t,u,v):
 up=np.stack([np.ones_like(u),u,u*u],-1);vp=np.stack([np.ones_like(v),v,v*v],-1)
 prof=np.stack([up[...,i]*vp[...,j]*(1 if i==j else 1)+ (up[...,j]*vp[...,i] if i!=j else 0)for i,j in ij],-1)
 return (gegen(t)[...,None]*weights[:,None]*prof[...,None,:]).reshape(-1,N*6)
def trirows(x):
 a=1+6/13*x[:,0];b=a+6/13*(1-x[:,0])*x[:,1];c=b+6/13*(1-x[:,0])*(1-x[:,1])*x[:,2];u=(1-x[:,0])**4;w=u*(1-x[:,1])**4;z=w*(1-x[:,2])**4
 return bilrows((a+b-c)/(2*np.sqrt(a*b)),u,w)+(1-x[:,2,None])**p*bilrows((a+c-b)/(2*np.sqrt(a*c)),u,z)+((1-x[:,1,None])*(1-x[:,2,None]))**p*bilrows((b+c-a)/(2*np.sqrt(b*c)),w,z)
def edge_points(vec):
 points=[]
 for axis in range(3):
  for fixed in itertools.product([0.,1.],repeat=2):
   grid=np.linspace(0,1,257);pts=np.tile(np.insert(fixed,axis,0.),(len(grid),1));pts[:,axis]=grid;vals=trirows(pts)@vec;points.extend(pts[[0,-1]])
   for j in range(1,len(grid)-1):
    if vals[j]>=vals[j-1] and vals[j]>=vals[j+1]:
     def fun(x):
      pt=pts[j].copy();pt[axis]=x;return -(trirows(pt[None,:])@vec).item()
     opt=minimize_scalar(fun,bounds=(grid[j-1],grid[j+1]),method='bounded',options={'xatol':1e-13});pt=pts[j].copy();pt[axis]=opt.x;points.append(pt)
   points.extend(pts[np.argsort(vals)[-10:]])
 return np.array(points)
def scale(rows):return rows/np.maximum(np.max(abs(rows),axis=1,keepdims=True),1e-200)
rng=np.random.default_rng(seed);xs=np.r_[rng.random((2500,3)),np.array(np.meshgrid(*[[0,.001,.01,.05,.1,.2,.3,.5,.7,.9,.99,1]]*3)).reshape(3,-1).T]
cutfile=D/f'guarded_cuts_p{p}.json'
if cutfile.exists():xs=np.r_[xs,np.array(json.loads(cutfile.read_text()))]
xs=np.r_[xs,[[0,0,.625]]]
T=np.array([[1,0,0],[1,.5,0],[1,1,1.]])
capfile=D/f'cap_cuts_p{p}.json'
t=np.r_[np.linspace(-1,.5,501),json.loads(capfile.read_text())if capfile.exists()else []];gg=gegen(t)*weights;caps=[]
for a,b in ij:
 prof=np.array([T[a,i]*T[b,j]+(T[a,j]*T[b,i]if i!=j else 0)for i,j in ij]);caps.append((gg[:,:,None]*prof).reshape(-1,N*6))
cap=scale(np.vstack(caps));oldcaps=json.load(open(D/'base/scalar/caps/result.json'))
raw=json.load(open(D/'base/scalar/repaired.json'));P=arb_poly([arb(fmpq(x))for x in raw['P']]);points=np.linspace(22,22*19/13,1401);fv=np.array([float((P(arb(float(x)))*arb(-float(x)).exp()).mid())for x in points]);h=(19/13-points/22)/(6/13);H=np.stack([h**p,h**(p+4),h**(p+8)],-1);diag=np.tile(np.stack([H[:,i]*H[:,j]*(1 if i==j else 2)for i,j in ij],-1),(1,N))*np.repeat(weights,6)
allow=[]
for x in points:
 bins=[r for r in oldcaps['rows']if float(F(r['l']))-1e-10<=x<=float(F(r['r']))+1e-10];allow.append(min(float(F(r['height']))for r in bins))
allow=np.array(allow);X=cp.Variable((4,4),symmetric=True);Tleg=np.array([[1.,0,0],[-1,2,0],[1,-6,6]]); Cs=[X[1:,1:]]+[cp.Variable((3,3),symmetric=True)for _ in range(N-1)];Bs=[Tleg.T@c@Tleg for c in Cs];vec=cp.hstack([cp.hstack([b[i,j]for i,j in ij])for b in Bs]);z=X[0,1:]@Tleg;s=X[0,0]
degree=71;alpha=n/2-1;norm=np.exp(gammaln(np.arange(degree+1)+n/2)-gammaln(n/2)-gammaln(np.arange(degree+1)+1));sign=(-1.)**np.arange(degree+1)
def lag(x):return np.stack([eval_genlaguerre(j,alpha,2*np.asarray(x))/norm[j]for j in range(degree+1)],-1)
delta=cp.Variable(degree+1)
Qpoly=arb_poly([arb(fmpq(x))for x in raw['Q']]);rootpoints=[float(F(f['mid']))for f in raw['factors']if f['type']=='positive_double']
fgrid=np.unique(np.r_[np.linspace(0,500,2001),np.geomspace(500,2500,151),rootpoints,[x+e for x in rootpoints for e in [-.01,-.0001,.0001,.01]]]);qval=np.array([float(Qpoly(arb(float(x))).mid())for x in fgrid]);fmat=lag(fgrid)*sign;fscale=np.maximum(np.maximum(np.max(abs(fmat),axis=1),abs(qval)),1e-50)
physical=np.r_[np.linspace(22*19/13,180,1601),np.linspace(180,500,501),np.geomspace(500,2500,151)];pval=np.array([float(P(arb(float(x))).mid())for x in physical]);pmat=lag(physical);pallow=[]
for xx in physical:
 bins=[r for r in oldcaps['rows']if float(F(r['l']))-1e-10<=xx<=float(F(r['r']))+1e-10]
 pallow.append(min(float(F(r['height']))for r in bins)*np.exp(xx)if bins else 0.)
pallow=np.array(pallow);pscale=np.maximum(np.maximum(np.max(abs(pmat),axis=1),abs(pval)),1e-50);pairmat=q*np.exp(-points[:,None])*lag(points)
for step in range(8):
 tr=scale(trirows(xs));gg=gegen(t)*weights;caps=[]
 for a,b in ij:
  prof=np.array([T[a,i]*T[b,j]+(T[a,j]*T[b,i]if i!=j else 0)for i,j in ij]);caps.append((gg[:,:,None]*prof).reshape(-1,N*6))
 cap=scale(np.vstack(caps));radscale=np.maximum(np.maximum(abs(fv*q),np.max(abs(diag),axis=1)),np.max(2*abs(H),axis=1));radscale=np.maximum(radscale,np.max(abs(pairmat),axis=1));radscale=np.maximum(radscale,1e-12)
 cons=[X>>1e-7*np.eye(4)]+[b>>1e-7*np.eye(3)for b in Cs[1:]]+[tr@vec<=-2e-6,cap@vec<=-2e-6,(diag/radscale[:,None])@vec+(pairmat/radscale[:,None])@delta+2*(H/radscale[:,None])@z+(fv-allow)*q/radscale<=-1e-8]
 cons.extend([sign@delta==0,(fmat/fscale[:,None])@delta+qval/fscale>=0,(pmat/pscale[:,None])@delta+(pval-pallow)/pscale<=0])
 pr=cp.Problem(cp.Minimize(s+cp.sum(delta)),cons);start=time.time()
 try:pr.solve(solver='CLARABEL',max_iter=350,max_threads=2,tol_gap_abs=1e-9,tol_feas=1e-9,tol_gap_rel=1e-9)
 except Exception as e:
  print(repr(e),flush=True);break
 if X.value is None:print(pr.status,flush=True);break
 (D/f'guarded_raw_p{p}_step{step}.json').write_text(json.dumps(dict(status='NUMERICAL_UNSCREENED',s=float(s.value),z=z.value.tolist(),B=[b.value.tolist()for b in Bs]),indent=2));extreme=differential_evolution(lambda x:-(trirows(np.asarray(x)[None,:])@vec.value).item(),[(0,1)]*3,maxiter=150,popsize=15,tol=1e-9,seed=51);test=np.r_[rng.random((30000,3)),xs,extreme.x[None,:],edge_points(vec.value)];res=trirows(test)@vec.value;bad=np.argsort(res)[-30:];record=dict(status='NUMERICAL_ONLY_NOT_A_BOUND',radial_delta=delta.value.tolist(),objective_delta=float(s.value+sum(delta.value)),dimension=n,powers=[p,p+4,p+8],iteration=step,solver_status=pr.status,s=float(s.value),z=z.value.tolist(),B=[b.value.tolist()for b in Bs],q=q,co=co.tolist(),R='19/13',maximum_sampled_triangle=float(res.max()),minimum_eigenvalue=min(np.linalg.eigvalsh(X.value).min(),*[np.linalg.eigvalsh(b.value).min()for b in Bs]),seconds=time.time()-start,source_commit='34f59648d034e92a52768c95af96221a8f589ac9',old_s=float(F(v['s'])))
 (D/f'joint_p{p}_step{step}.json').write_text(json.dumps(record,indent=2));print({key:val for key,val in record.items()if key not in ['B','z','co']},flush=True)
 newcap,capmax=extrema(v['co'],[b.value.tolist()for b in Bs]);print('EXACT CAP EXTREMUM',capmax,flush=True);t=np.unique(np.r_[t,newcap]);capfile.write_text(json.dumps(t.tolist()))
 if res.max()<0 and capmax<0:break
 xs=np.r_[xs,test[bad]];(D/f'guarded_cuts_p{p}.json').write_text(json.dumps(xs.tolist()))
