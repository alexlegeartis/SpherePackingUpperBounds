"""Discover deep holes by LP cutting planes for the Voronoi cell.
The CVP oracle uses floating fplll; no certified density claim is made.
"""
from pathlib import Path
import argparse,json,time,warnings
import numpy as np
from scipy.optimize import linprog,OptimizeWarning
from fpylll import Enumeration,GSO,IntegerMatrix
warnings.simplefilter('ignore',OptimizeWarning)
P=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--seed',type=int,default=13);ap.add_argument('--seconds',type=float,default=600);ap.add_argument('--starts',type=int,default=30);ap.add_argument('--batch',type=int,default=1);args=ap.parse_args()
j=np.load(P/'contacts_reduced.npz');Q=j['Q'];n=36;V=j['V'];V=np.unique(np.concatenate([V,-V,np.eye(n,dtype=np.int64),-np.eye(n,dtype=np.int64)]),axis=0);A=V@Q;b=np.einsum('ij,ij->i',A,V)/2
M=GSO.Mat(IntegerMatrix.from_matrix(Q.tolist()),gram=True,update=True);mu=np.array([[M.get_mu(j,i) if j>i else int(j==i) for i in range(n)]for j in range(n)])
R=np.linalg.cholesky(Q);rng=np.random.default_rng(args.seed);start=time.time();best=0;history=[];cuts=0
for restart in range(args.starts):
 if time.time()-start>args.seconds:break
 objective=rng.normal(size=n)@R.T
 for outer in range(15):
  for step in range(1000):
   if time.time()-start>args.seconds:break
   r=linprog(-objective,A_ub=A,b_ub=b,bounds=[(None,None)]*n,method='highs',options={'threads':1,'time_limit':30})
   if not r.success:print('LP fail',r.message,flush=True);break
   x=r.x;qx=x@Q@x
   E=Enumeration(M,nr_solutions=args.batch);res=E.enumerate(0,n,max(qx+1,180),0,target=(x@mu).tolist());z=np.rint(res[0][1]).astype(np.int64);d=(x-z)@Q@(x-z)
   if d>=qx-1e-6:
    if qx>best+1e-5:
     best=qx;record={'distance_squared':qx/20,'shift_coefficients':x.tolist(),'nearest_oracle_squared':d/20,'cuts':cuts,'restart':restart,'outer':outer,'seconds':time.time()-start,'scope':'numerical LP and floating CVP discovery, not certificate'};history.append(record);(P/f'voronoi_seed{args.seed}.json').write_text(json.dumps({'best':record,'history':history},indent=2));print(json.dumps(record),flush=True)
    break
   if not np.any(z):print('NUMERICAL ORACLE BUG',flush=True);raise RuntimeError
   zs=np.array([np.rint(rr[1]).astype(np.int64) for rr in res if rr[0]<qx-1e-6]); az=zs@Q;rhs=np.einsum('ij,ij->i',zs,az)/2;A=np.vstack([A,az,-az]);b=np.concatenate([b,rhs,rhs]);cuts+=len(zs)
   if args.batch>1 or cuts%20==0:print('cuts',cuts,'qx',qx/20,'distance',d/20,'time',time.time()-start,flush=True)
  else:break
  if time.time()-start>args.seconds or not r.success:break
  newobj=x@Q
  if np.linalg.norm(newobj/np.linalg.norm(newobj)-objective/np.linalg.norm(objective))<1e-6:break
  objective=newobj
 print('RESTART',restart,'best',best/20,'cuts',cuts,'time',time.time()-start,flush=True)
print('DONE best',best/20,'cuts',cuts,'time',time.time()-start,flush=True)
