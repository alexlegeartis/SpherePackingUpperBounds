"""Numerical discovery of a joint 37D metric from a CT36 deep-hole seed.
No result from this script is an exact shortest-vector certificate.
"""
import argparse,json,time,warnings,math
from pathlib import Path
import numpy as np
from scipy.linalg import qr
from scipy.optimize import linprog,OptimizeWarning
from fpylll import Enumeration,GSO,IntegerMatrix
import sympy as sp
warnings.simplefilter('ignore',OptimizeWarning)
P=Path(__file__).resolve().parent;ROOT=P.parents[2]
ap=argparse.ArgumentParser();ap.add_argument('--seconds',type=float,default=420);ap.add_argument('--iterations',type=int,default=30);ap.add_argument('--step',type=float,default=.015);args=ap.parse_args();start=time.time()
Q=np.loadtxt(ROOT/'constructions/dim36/Q_reduced.txt',dtype=np.int64)
j=max([json.loads(p.read_text())['best'] for p in P.glob('voronoi_seed*.json')],key=lambda r:r['distance_squared']);x=np.array(j['shift_coefficients']);r=x@Q@x
M=GSO.Mat(IntegerMatrix.from_matrix(Q.tolist()),gram=True,update=True);mu=np.array([[M.get_mu(j,i) if j>i else int(j==i) for i in range(36)]for j in range(36)])
E=Enumeration(M,nr_solutions=500);near=E.enumerate(0,36,r+1e-5,0,target=(x@mu).tolist());Z=np.array([np.rint(a[1]).astype(np.int64) for a in near]);Z=Z[np.any(Z,axis=1)]
_,R,piv=qr(Z.T,pivoting=True,mode='economic');rank=np.linalg.matrix_rank(Z);assert rank==36,rank;T=Z[piv[:36]];A=sp.Matrix(T@Q);b=sp.Matrix([int(z@Q@z)//2 for z in T]);xr=A.inv()*b
xf=np.array(xr,dtype=float).reshape(-1);assert np.max(abs(xf-x))<1e-6
x=xf;cross=Q@x/20;H=np.block([[Q/20,cross[:,None]],[cross[None,:],np.array([[8.]])]])
H0=H.copy();n=37;ii,jj=np.triu_indices(n);factors=np.where(ii==jj,1.,2.);obj=np.where(ii==jj,1.,0.)
raw=np.loadtxt(ROOT/'research/density_breakthroughs_20261113/contact_walk/contacts_reduced.txt',dtype=np.int64);V=np.pad(raw[:,1:],((0,0),(0,1)));W=np.column_stack([-Z,np.ones(len(Z),dtype=np.int64)]);W=np.vstack([W,np.eye(n,dtype=np.int64)[-1]])
rng=np.random.default_rng(37);select=rng.choice(len(V),min(4000,len(V)),replace=False);cuts=np.vstack([V[select],W]);cutkeys={tuple(v) for v in cuts};best=0.;trace=[]
def density(H):return math.exp(n/2*math.log(2)-np.linalg.slogdet(H)[1]/2)
def save(H,N,m,it,status):
 global best
 d=density(H)
 if d>best:
  best=d;np.savez(P/'joint37_best.npz',H=H,integer_Gram=N,minimum_numerator=m,scale=np.array(100000000))
  row={'iteration':it,'center_density':d,'record_center_density':4*math.sqrt(2),'ratio_to_record':d/(4*math.sqrt(2)),'minimum_numerator':int(m),'scale':100000000,'cuts':len(cuts),'seconds':time.time()-start,'status':status,'scope':'floating CVP/SVP discovery, not a minimum certificate'};trace.append(row);(P/'joint37_results.json').write_text(json.dumps({'best':row,'trace':trace,'seed_hole':j,'seed_shift_rational':[str(t) for t in xr]},indent=2));print(json.dumps(row),flush=True)
print('SEED',density(H),'hole',x@Q@x/20,'nearest',len(Z)+1,'contacts',len(cuts),flush=True)
save(H,np.rint(H*1e8).astype(np.int64),800000000,-1,'unverified seed')
for it in range(args.iterations):
 if time.time()-start>args.seconds:break
 L=np.linalg.cholesky(H);accepted=False
 for trial in range(60):
  if time.time()-start>args.seconds:break
  coords=cuts@L;norms=np.einsum('bi,bi->b',coords,coords);C=coords[:,ii]*coords[:,jj]*factors
  sol=linprog(obj,A_ub=-C,b_ub=norms-8.,bounds=[(-args.step,args.step)]*len(ii),method='highs',options={'threads':1,'time_limit':40,'dual_feasibility_tolerance':1e-8,'primal_feasibility_tolerance':1e-8})
  if not sol.success:print('LPFAIL',sol.message,flush=True);break
  print('LP',it,trial,'trace',sol.fun,'cuts',len(cuts),'time',time.time()-start,flush=True)
  if sol.fun>-1e-7:break
  D=np.zeros((n,n));D[ii,jj]=sol.x;D[jj,ii]=sol.x
  alpha=min(1.,.5/max(1e-10,-np.linalg.eigvalsh(D)[0]));candidate=L@(np.eye(n)+alpha*D)@L.T
  vals=np.einsum('bi,ij,bj->b',V,candidate,V);bad=np.flatnonzero(vals<8.-1e-7)
  if len(bad):
   use=bad[np.argsort(vals[bad])[:200]];new=[V[t] for t in use if tuple(V[t]) not in cutkeys]
   if not new:print('BAD KNOWN CONTACT',float(vals.min()),flush=True);break
   cuts=np.vstack([cuts,new]);cutkeys.update(map(tuple,new));print('CORECUT',len(new),len(bad),'min',float(vals.min()),flush=True);continue
  N=np.rint(candidate*1e8).astype(np.int64);B=IntegerMatrix.from_matrix(N.tolist());MG=GSO.Mat(B,gram=True,update=True);EE=Enumeration(MG);out=EE.enumerate(0,n,800000100,0);z=np.rint(out[0][1]).astype(np.int64);m=int(z@N@z)
  if m<799999900:
   if tuple(z) in cutkeys or tuple(-z) in cutkeys:print('BAD KNOWN SVP',m,flush=True);break
   cuts=np.vstack([cuts,z]);cutkeys.add(tuple(z));print('SVPCUT',m,'nodes',EE.get_nodes(),flush=True);continue
  H=N.astype(float)*(8./m);accepted=True;save(H,N,m,it,'floating SVP completed');break
 if not accepted:break
print('DONE','best',best,'seconds',time.time()-start,flush=True)
