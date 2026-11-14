"""Numerical deep-hole search for a prospective laminated 37D packing.

Every reported distance is an UPPER bound from an explicit nearest-point
witness. No successful search or CVP result is a minimum certificate.
"""
from pathlib import Path
import json,time
import numpy as np
from scipy.optimize import linprog
from fpylll import IntegerMatrix,LLL,CVP

P=Path(__file__).resolve().parent
Q=np.load(P/'new36.npz')['Q'];G=Q/20;n=36
V=np.loadtxt(P/'new36_vectors.txt',dtype=np.int64)
assert np.all(np.einsum('vi,ij,vj->v',V,Q,V)==160)
A=V@G;h=np.full(len(V),4.)
B=np.rint(np.linalg.cholesky(G)*1e7).astype(np.int64)
b=IntegerMatrix.from_matrix(B.tolist());U=IntegerMatrix.identity(n);LLL.reduction(b,U=U)
Br=np.array([[int(b[i,j])for j in range(n)]for i in range(n)],dtype=np.int64)
Ui=np.array([[int(U[i,j])for j in range(n)]for i in range(n)],dtype=np.int64);Bi=np.linalg.inv(Br)
def cvp(x):
 v=np.array(CVP.closest_vector(b,tuple(np.rint(x@B).astype(np.int64)),method='fast'),dtype=np.int64)
 return np.rint(v@Bi).astype(np.int64)@Ui

rng=np.random.default_rng(20261114);seeds=np.load(P/'cvp37_seeds.npz');ids=np.argsort(seeds['norms'])[::-1][:12]
rows=[];start=time.time()
for trial,idx in enumerate(ids):
 x=seeds['X'][idx]-seeds['Z'][idx]
 added=[];added_h=[];solves=0
 active=rng.choice(len(V),2000,replace=False)
 for step in range(5):
  c=-(x@G)
  for cut in range(50):
   res=linprog(c,A_ub=np.vstack([A[active],np.array(added).reshape(-1,n)]),b_ub=np.r_[h[active],added_h],bounds=[(None,None)]*n,method='highs',options={'time_limit':15.})
   solves+=1
   if not res.success:break
   y=res.x
   violations=A@y-h
   bad=np.flatnonzero(violations>1e-7)
   if len(bad):
    active=np.unique(np.r_[active,bad[np.argsort(violations[bad])[-250:]]])
    continue
   z=cvp(y);d=float((y-z)@G@(y-z));r=float(y@G@y)
   if r-d<1e-7:break
   added.append(z@G);added_h.append(float(z@G@z)/2)
  if not res.success:break
  delta=float(np.linalg.norm(y-x));x=y
  if delta<1e-7:break
 z=cvp(x);dist=float((x-z)@G@(x-z))
 rows.append({'trial':trial,'point':x.tolist(),'nearest_witness':z.tolist(),'squared_distance_upper':dist,'cuts':len(added),'lp_solves':solves,'seconds':time.time()-start})
 (P/'holes37_results.json').write_text(json.dumps(rows,indent=2)+'\n')
 print('trial',trial,'distance upper',dist,'cuts',len(added),'solves',solves,'seconds',time.time()-start,flush=True)
