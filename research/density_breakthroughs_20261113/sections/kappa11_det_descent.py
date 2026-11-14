"""Integral Gram discovery by log-determinant linearization and SVP cuts.

This only proposes small normal lattices; no ambient embedding is asserted.
The floating enumerator is a discovery oracle, not a minimum certificate.
"""
from pathlib import Path
import numpy as np, json, subprocess, time, sys
from scipy.optimize import milp, Bounds, LinearConstraint
import sympy as s
P=Path(__file__).resolve().parent;ROOT=P.parents[2];n=11;i,j=np.triu_indices(n,1)
raw=json.load(open(P/'kappa11_integer_neighbors.json'))
H=np.array(min(raw,key=lambda r:r['determinant'])['gram'],dtype=np.int64)
old=np.loadtxt(P/'kappa11_vectors.txt',dtype=np.int64)[:,1:]
vectors={tuple(v) for v in old}
out=[];start=time.time();seconds=float(sys.argv[1])if len(sys.argv)>1 else 120.
rng=np.random.default_rng(11)
def shell(Q,threshold):
 q=P/'integer_descent_Q.txt';r=P/'integer_descent_R.txt';v=P/'integer_descent_vectors.txt'
 np.savetxt(q,Q,fmt='%d');np.savetxt(r,np.linalg.cholesky(Q).T,fmt='%.17g')
 run=subprocess.run([str(ROOT.parent/'shell_discovery'),str(n),str(q),str(r),str(threshold),'1000000',str(v),'0'],text=True,capture_output=True)
 stats=json.loads(run.stdout);assert not stats['capped']
 if not v.stat().st_size:return np.empty((0,n),dtype=np.int64)
 rows=np.loadtxt(v,dtype=np.int64,ndmin=2)
 assert np.all(np.einsum('bi,ij,bj->b',rows[:,1:],Q,rows[:,1:])==rows[:,0])
 return rows[:,1:]
for v in shell(H,8):vectors.add(tuple(v))
best=int(s.Matrix(H).det());bestH=H.copy();iterations=0;unchanged=0
while time.time()-start<seconds:
 V=np.array(list(vectors),dtype=np.int64);A=2*V[:,i]*V[:,j];b=6-6*np.square(V).sum(1)
 inv=np.linalg.inv(H.astype(float));obj=2*inv[i,j]
 if iterations%4==3:obj+=rng.normal(size=len(i))*.005
 r=milp(obj,integrality=np.ones(len(i)),bounds=Bounds(np.full(len(i),-3),np.full(len(i),3)),constraints=LinearConstraint(A,b,np.inf),options={'time_limit':min(15,max(1,seconds-(time.time()-start))),'mip_rel_gap':0})
 if r.x is None:break
 C=6*np.eye(n,dtype=np.int64);x=np.rint(r.x).astype(np.int64);C[i,j]=x;C[j,i]=x
 eig,Q=np.linalg.eigh(C.astype(float));row={'iteration':iterations,'constraints':len(V),'solver':r.message,'elapsed':time.time()-start,'minimum_eigenvalue':float(eig[0])}
 if eig[0]<=1e-8:
  # Any exact nonpositive integer direction supplies a necessary lattice cut.
  added=0
  for scale in range(2,500):
   z=np.rint(scale*Q[:,0]).astype(np.int64)
   if np.any(z) and z@C@z<6:
    vectors.add(tuple(z));added+=1
    if added>=4:break
  assert added
  row['indefinite_cuts']=added
 else:
  W=shell(C,8);norms=np.einsum('bi,ij,bj->b',W,C,W)
  for v in W:vectors.add(tuple(v))
  row['short_vector_count']=int(np.sum(norms<6));row['determinant']=int(s.Matrix(C).det())
  if row['short_vector_count']==0:
   unchanged=unchanged+1 if np.array_equal(H,C) else 0
   H=C
   if row['determinant']<best:
    best=row['determinant'];bestH=C.copy()
    np.savez(P/'integer_rank11_best.npz',G=C)
    row['best']=best
  else:
   row['short_vector']=W[np.argmin(norms)].tolist()
 out.append(row);print(json.dumps(row),flush=True)
 (P/'kappa11_det_descent.json').write_text(json.dumps({'best':best,'gram':bestH.tolist(),'rows':out},indent=2))
 iterations+=1
 if unchanged>=16:break
print('DONE',best,'iterations',iterations,'seconds',time.time()-start,flush=True)
