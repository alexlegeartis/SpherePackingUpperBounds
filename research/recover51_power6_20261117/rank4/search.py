"""Four-profile compact PSD kernel reconstruction. Numerical research, NOT a bound."""
import argparse,json,sys,time,math,itertools,hashlib
from pathlib import Path
from fractions import Fraction as F
import numpy as np
import cvxpy as cp
from scipy.special import eval_gegenbauer
from scipy.optimize import differential_evolution,minimize_scalar
from flint import arb,arb_poly,ctx,fmpq as Q,fmpq_poly as Poly
sys.set_int_max_str_digits(0);ctx.prec=512
ap=argparse.ArgumentParser();ap.add_argument('--power',type=int,default=4);ap.add_argument('--spacing',type=int,default=4);ap.add_argument('--rounds',type=int,default=8);a=ap.parse_args()
ROOT=Path(__file__).parent;D=ROOT/f'p{a.power}_d{a.spacing}';D.mkdir(exist_ok=True);BASE=ROOT.parent/'base';p=a.power;spacing=a.spacing;r=4;m=r*(r+1)//2;powers=[p+spacing*i for i in range(r)];seed=5104
v=json.load(open(BASE/'geometry/dim51/frozen.json'));n=51;q=float(F(v['q']));co=np.array([float(F(x))for x in v['co']]);weights=co/q;N=len(co);kk=np.arange(N);den=eval_gegenbauer(kk,(n-2)/2,1);ij=list(itertools.combinations_with_replacement(range(r),2));rng=np.random.default_rng(seed)
T=np.array([[math.comb(j,i)/math.comb(r-1,i) if i<=j else 0 for i in range(r)]for j in range(r)])
L=np.array([[1.,0,0,0],[-1,2,0,0],[1,-6,6,0],[-1,12,-30,20.]])
def gegen(t):return eval_gegenbauer(kk,(n-2)/2,np.asarray(t)[...,None])/den

def bilrows(t,u,v):
 U=np.stack([u**j for j in range(r)],-1);V=np.stack([v**j for j in range(r)],-1);prod=np.stack([U[...,i]*V[...,j]+(U[...,j]*V[...,i]if i!=j else 0)for i,j in ij],-1)
 return (gegen(t)[...,None]*weights[:,None]*prod[...,None,:]).reshape(-1,N*m)
def trirows(x):
 x=np.atleast_2d(x);aa=1+6/13*x[:,0];bb=aa+6/13*(1-x[:,0])*x[:,1];cc=bb+6/13*(1-x[:,0])*(1-x[:,1])*x[:,2];u=(1-x[:,0])**spacing;w=u*(1-x[:,1])**spacing;z=w*(1-x[:,2])**spacing
 return bilrows((aa+bb-cc)/(2*np.sqrt(aa*bb)),u,w)+(1-x[:,2,None])**p*bilrows((aa+cc-bb)/(2*np.sqrt(aa*cc)),u,z)+((1-x[:,1,None])*(1-x[:,2,None]))**p*bilrows((bb+cc-aa)/(2*np.sqrt(bb*cc)),w,z)
def scale(x):return x/np.maximum(np.max(abs(x),axis=1,keepdims=True),1e-200)
def boundary_grid(size=17):
 pts=[]
 for axis in range(3):
  for endpoint in [0.,1.]:
   grid=np.array(np.meshgrid(np.linspace(0,1,size),np.linspace(0,1,size))).reshape(2,-1).T;pts.append(np.insert(grid,axis,endpoint,axis=1))
 return np.concatenate(pts)
def separate(vec):
 out=[]
 for axis in range(3):
  for fixed in itertools.product([0.,1.],repeat=2):
   grid=np.linspace(0,1,513);pts=np.tile(np.insert(fixed,axis,0.),(len(grid),1));pts[:,axis]=grid;vals=trirows(pts)@vec;out.extend(pts[np.argsort(vals)[-10:]])
   for j in range(1,len(grid)-1):
    if vals[j]>=vals[j-1] and vals[j]>=vals[j+1]:
     def f(t):
      pt=pts[j].copy();pt[axis]=t;return -(trirows(pt)@vec).item()
     opt=minimize_scalar(f,bounds=(grid[j-1],grid[j+1]),method='bounded',options={'xatol':1e-13});pt=pts[j].copy();pt[axis]=opt.x;out.append(pt)
 for axis in range(3):
  for endpoint in [0.,1.]:
   def f(z):return -(trirows(np.insert(z,axis,endpoint))@vec).item()
   opt=differential_evolution(f,[(0,1)]*2,maxiter=80,popsize=10,seed=seed+axis,tol=1e-9);out.append(np.insert(opt.x,axis,endpoint))
 opt=differential_evolution(lambda x:-(trirows(x)@vec).item(),[(0,1)]*3,maxiter=150,popsize=12,seed=seed,tol=1e-9);out.append(opt.x)
 return np.array(out)
def cap_extrema(B):
 coq=list(map(Q,v['co']));B=[[[Q(str(F(str(x))))for x in row]for row in b]for b in B];x=Poly([0,1]);G=[Poly([1]),x]
 for k in range(2,N):G.append(((2*k+n-4)*x*G[-1]-(k-1)*G[-2])/(k+n-3))
 Tq=[[Q(math.comb(j,i),math.comb(r-1,i))if i<=j else Q(0)for i in range(r)]for j in range(r)];pts=[];mx=arb(-10**100)
 for aa,bb in ij:
  poly=sum((coq[k]*sum((Tq[aa][i]*B[k][i][j]*Tq[bb][j]for i in range(r)for j in range(r)),Q(0))*G[k]for k in range(N)),Poly([0]));ap=arb_poly([arb(x)for x in poly]);roots=poly.derivative().complex_roots();assert sum(m for z,m in roots)==poly.degree()-1;mx=max(mx,arb(poly(-1)).upper(),arb(poly(Q(1,2))).upper())
  for z,mult in roots:
   if z.imag.contains(0) and not z.real < -1 and not z.real > arb(1)/2:pts.append(float(z.real.mid()));mx=max(mx,ap(z.real).upper())
 return pts,float(mx)
xs=np.r_[rng.random((1800,3)),boundary_grid(),np.array(np.meshgrid(*[[0,.01,.1,.25,.5,.625,.75,.9,.99,1]]*3)).reshape(3,-1).T]
cutfile=D/'tri_cuts.json';capfile=D/'cap_cuts.json'
if cutfile.exists():xs=np.r_[xs,np.array(json.load(open(cutfile)))];xs=np.unique(xs,axis=0)
ts=np.r_[np.linspace(-1,.5,601),json.load(open(capfile))if capfile.exists()else []]
oldcaps=json.load(open(BASE/'scalar/caps/result.json'));raw=json.load(open(BASE/'scalar/repaired.json'));P=arb_poly([arb(Q(x))for x in raw['P']]);points=np.linspace(22,22*19/13,1401);fv=np.array([float((P(arb(float(x)))*arb(-float(x)).exp()).mid())for x in points]);h=(19/13-points/22)/(6/13);H=np.stack([h**z for z in powers],-1);diag=np.tile(np.stack([H[:,i]*H[:,j]*(1 if i==j else 2)for i,j in ij],-1),(1,N))*np.repeat(weights,m)
allow=[]
for x in points:
 rows=[r for r in oldcaps['rows']if float(F(r['l']))-1e-10<=x<=float(F(r['r']))+1e-10];allow.append(min(float(F(r['height']))for r in rows))
allow=np.array(allow);X=cp.Variable((r+1,r+1),symmetric=True);Cs=[X[1:,1:]]+[cp.Variable((r,r),symmetric=True)for _ in range(N-1)];Bs=[L.T@c@L for c in Cs];vec=cp.hstack([cp.hstack([b[i,j]for i,j in ij])for b in Bs]);z=X[0,1:]@L;s=X[0,0]
radscale=np.maximum(np.maximum(abs(fv*q),np.max(abs(diag),axis=1)),np.max(2*abs(H),axis=1));radscale=np.maximum(radscale,1e-12)
for step in range(a.rounds):
 start=time.time();tr=scale(trirows(xs));gg=gegen(ts)*weights;caps=[]
 for aa,bb in ij:
  prof=np.array([T[aa,i]*T[bb,j]+(T[aa,j]*T[bb,i]if i!=j else 0)for i,j in ij]);caps.append((gg[:,:,None]*prof).reshape(-1,N*m))
 cap=scale(np.vstack(caps));cons=[X>>1e-7*np.eye(r+1)]+[c>>1e-7*np.eye(r)for c in Cs[1:]]+[tr@vec<=-2e-6,cap@vec<=-2e-6,(diag/radscale[:,None])@vec+2*(H/radscale[:,None])@z+(fv-allow)*q/radscale<=-1e-8]
 print('SOLVING',powers,step,'triangles',len(xs),'cap_rows',len(cap),flush=True);pr=cp.Problem(cp.Minimize(s),cons)
 try:pr.solve(solver='CLARABEL',max_iter=350,max_threads=2,tol_gap_abs=1e-9,tol_feas=1e-9,tol_gap_rel=1e-9)
 except Exception as e:(D/'failure.json').write_text(json.dumps({'status':'SOLVER_FAILURE','error':repr(e),'iteration':step}));print(repr(e),flush=True);break
 if X.value is None:print(pr.status,flush=True);break
 record=dict(status='NUMERICAL_ONLY_NOT_A_BOUND',dimension=n,powers=powers,iteration=step,solver_status=pr.status,s=float(s.value),z=z.value.tolist(),B=[b.value.tolist()for b in Bs],q=q,co=co.tolist(),R='19/13',old_s=float(F(v['s'])),source_commit='34f59648d034e92a52768c95af96221a8f589ac9',source_archive_sha256='7d1b48627948429301d916f3ed04beb255d33c7cd325cb47a7ca81acbb7128e6',minimum_eigenvalue=float(min(np.linalg.eigvalsh(X.value).min(),*[np.linalg.eigvalsh(b.value).min()for b in Bs])),seconds=time.time()-start)
 (D/f'unscreened{step}.json').write_text(json.dumps(record,indent=2));test=np.r_[rng.random((20000,3)),xs,boundary_grid(41),separate(vec.value)];res=trirows(test)@vec.value;bad=np.argsort(res)[-50:];newcaps,capmax=cap_extrema(record['B']);record.update(maximum_sampled_triangle=float(res.max()),maximum_exact_cap=capmax,seconds=time.time()-start,worst_point=test[np.argmax(res)].tolist());(D/f'candidate{step}.json').write_text(json.dumps(record,indent=2));print({k:v for k,v in record.items()if k not in ['B','z','co']},flush=True)
 ts=np.unique(np.r_[ts,newcaps]);capfile.write_text(json.dumps(ts.tolist()));xs=np.unique(np.r_[xs,test[bad]],axis=0);cutfile.write_text(json.dumps(xs.tolist()))
 if res.max()<0 and capmax<0 and record['minimum_eigenvalue']>0:print('ALL DISCOVERY SCREENS PASS; FULL PROOFS AND INDEPENDENT AUDIT STILL REQUIRED',flush=True);break
