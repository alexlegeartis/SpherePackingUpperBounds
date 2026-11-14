"""Exploratory deep-hole search in seven determinant2^16 32D lattices.

CVP and LP results are floating discovery output, never minimum certificates.
Any hole farther than 4 from its parent gives a record candidate in dimension33
via a two-periodic layer construction (without requiring a torsion shift).
"""
from pathlib import Path
import argparse,json,time,subprocess,warnings,os
import numpy as np
import sympy as sp
from scipy.optimize import linprog,OptimizeWarning
from fpylll import Enumeration,GSO,IntegerMatrix

warnings.simplefilter('ignore',OptimizeWarning)
P=Path(__file__).resolve().parent
ROOT=P.parents[2]
ap=argparse.ArgumentParser();ap.add_argument('index',type=int);ap.add_argument('--seconds',type=float,default=90);ap.add_argument('--seed',type=int,default=1472);args=ap.parse_args()
i=args.index;n=32;suffix='_bkz' if (P/f'Q32_{i}_bkz.txt').exists() else ''
G=np.loadtxt(P/f'Q32_{i}{suffix}.txt',dtype=np.int64)
pool=P/f'minvecs_{i}{suffix}.txt'
if not pool.exists():
 s=subprocess.run([str(ROOT/'research/density_breakthroughs_20261113/bin/shell_discovery'),'32',str(P/f'Q32_{i}{suffix}.txt'),str(P/f'Q32_{i}{suffix}_R.txt'),'6','100000000',str(pool),'0'],capture_output=True,text=True,check=True)
 (P/f'minvecs_{i}.json').write_text(s.stdout)
V=np.loadtxt(pool,dtype=np.int64)[:,1:]
assert np.all(np.einsum('bi,ij,bj->b',V,G,V)==6)
M=GSO.Mat(IntegerMatrix.from_matrix(G.tolist()),gram=True,update=True)
mu=np.array([[M.get_mu(j,k) if j>k else int(j==k) for k in range(n)]for j in range(n)])
rng=np.random.default_rng(args.seed+i);rng.shuffle(V)
A0=V[:2500]@G;b0=np.full(len(A0),3.)
A=np.vstack([A0,-A0,np.eye(n)@G,-np.eye(n)@G]);b=np.r_[b0,b0,np.diag(G)/2,np.diag(G)/2]
R=np.linalg.cholesky(G);start=time.time();best=0.;cuts=0;calls=0;history=[];restarts=0
hole_seed=None
sf=P/f'half_holes_{i}_bkz.json'
if sf.exists():
 sj=json.loads(sf.read_text());hole_seed=np.array(sj['witness'],dtype=float)/2
else:
 sf=P/f'half_holes_{i}.json'
 if sf.exists():
  sj=json.loads(sf.read_text())
  if sj['best_found_norm']==16:
   hole_seed=np.array(sj['witness'],dtype=float)/2
   if suffix:
    T=np.load(P/f'Q32_{i}_bkz.npz')['T'];hole_seed=hole_seed@np.array(sp.Matrix(T).inv(),dtype=float)
if hole_seed is not None:
 ee=Enumeration(M).enumerate(0,n,6.,0,target=(hole_seed@mu).tolist())
 if abs(ee[0][0]-4)>1e-7:hole_seed=None
while time.time()-start<args.seconds:
 restarts+=1;objective=rng.normal(size=n)@R.T
 if restarts%3==0 and hole_seed is not None:objective=hole_seed@G+rng.normal(size=n)@R.T*rng.choice([.1,.3,.7,1.])
 elif restarts%3==0 and history:objective=np.array(history[-1]['x'])@G+rng.normal(size=n)@R.T*.2
 for outer in range(20):
  if time.time()-start>=args.seconds:break
  for step in range(500):
   if time.time()-start>=args.seconds:break
   r=linprog(-objective,A_ub=A,b_ub=b,bounds=[(None,None)]*n,method='highs',options={'threads':1,'time_limit':10})
   if not r.success:break
   x=r.x;qx=x@G@x;E=Enumeration(M,nr_solutions=12)
   res=E.enumerate(0,n,float(max(qx+1,6.)),0,target=(x@mu).tolist());calls+=1
   z=np.rint(res[0][1]).astype(np.int64);d=(x-z)@G@(x-z)
   if d>=qx-1e-7:
    if qx>best+1e-9:
     best=qx;record={'distance_squared':float(qx),'nearest_distance':float(d),'x':x.tolist(),'calls':calls,'cuts':cuts,'seconds':time.time()-start};history.append(record)
     print(json.dumps({'best':best,'calls':calls,'cuts':cuts,'seconds':time.time()-start}),flush=True)
     (P/f'voronoi32_{i}.json').write_text(json.dumps({'index':i,'history':history,'calls':calls,'cuts':cuts,'restarts':restarts,'scope':'Floating LP/CVP discovery; no rigorous minimum claim.'},indent=2))
    break
   zs=np.array([np.rint(rr[1]).astype(np.int64)for rr in res if rr[0]<qx-1e-7]);az=zs@G;rhs=np.einsum('ij,ij->i',zs,az)/2
   A=np.vstack([A,az,-az]);b=np.r_[b,rhs,rhs];cuts+=len(zs)
  if time.time()-start>=args.seconds or not r.success:break
  if best>4.000001:
   print('RECORD CANDIDATE',best,flush=True);break
  new=x@G
  if np.linalg.norm(new/np.linalg.norm(new)-objective/np.linalg.norm(objective))<1e-7:break
  objective=new
 if best>4.000001:break
report={'index':i,'basis_suffix':suffix,'best':best,'history':history,'calls':calls,'cuts':cuts,'restarts':restarts,'seconds':time.time()-start,'scope':'Floating LP/CVP discovery; no rigorous minimum claim.'}
(P/f'voronoi32_{i}.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items()if k!='history'}),flush=True)
