"""Search integral min-six forms near the scaled Kappa11 Gram.

Saved contact constraints are necessary only; a found candidate needs full SVP.
"""
import numpy as np,json,time
from pathlib import Path
from scipy.optimize import milp,Bounds,LinearConstraint
import sympy as s
P=Path(__file__).resolve().parent
G=np.array(json.load(open(P/'kappa11.json'))['gram']);target=1.5*G
raw=np.loadtxt(P/'kappa11_vectors.txt',dtype=np.int64);V=raw[:,1:];i,j=np.triu_indices(11,1)
A=2*V[:,i]*V[:,j];b=6-6*np.square(V).sum(1)
inv=np.linalg.inv(target);objective=2*inv[i,j]
out=[]
for width in [1,2,3]:
 lo=np.floor(target[i,j])-width;hi=np.ceil(target[i,j])+width
 r=milp(objective,integrality=np.ones(len(i)),bounds=Bounds(lo,hi),constraints=LinearConstraint(A,b,np.inf),options={'time_limit':40,'mip_rel_gap':0})
 row={'width':width,'status':r.message,'objective':float(r.fun)if r.fun is not None else None}
 if r.x is not None:
  H=6*np.eye(11,dtype=np.int64);x=np.rint(r.x).astype(np.int64);H[i,j]=x;H[j,i]=x
  assert np.all(np.einsum('bi,ij,bj->b',V,H,V)>=6)
  row.update(gram=H.tolist(),determinant=int(s.Matrix(H).det()),minimum_eigenvalue=float(np.linalg.eigvalsh(H)[0]),trace_ratio=float(np.trace(inv@H)))
 out.append(row);print(json.dumps(row),flush=True)
 (P/'kappa11_integer_neighbors.json').write_text(json.dumps(out,indent=2))
