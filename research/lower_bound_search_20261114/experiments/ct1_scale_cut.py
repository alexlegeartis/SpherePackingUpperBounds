"""Cutting-plane search in the other non-scalar CT/Turyn seed family."""
from pathlib import Path
import numpy as np,json,time
from scipy.spatial import HalfspaceIntersection
from fpylll import IntegerMatrix,LLL,BKZ
P=Path(__file__).resolve().parent;x=np.load(P/'CTturyn1_parts.npz');parts=x['parts'];den=int(x['denom']);Q0=x['Q'];V=np.eye(36,dtype=np.int64);co=np.array([np.einsum('vi,ij,vj->v',V,p,V)for p in parts]).T;co=np.unique(co,axis=0);rows=[];start=time.time()
for iteration in range(100):
 A=co/den;hs=np.r_[np.c_[-A,8*np.ones(len(A))],np.c_[-np.eye(3),np.full(3,0.0001)],np.c_[np.eye(3),np.full(3,-20.)]]
 poly=HalfspaceIntersection(hs,np.full(3,10.));vertices=poly.intersections;vertices=vertices[np.all(vertices>0,axis=1)];scale=min(vertices,key=np.prod);G=np.einsum('a,aij->ij',scale,parts)/den
 B=np.rint(np.linalg.cholesky(G)*1e7).astype(np.int64);b=IntegerMatrix.from_matrix(B.tolist());u=IntegerMatrix.identity(36);BKZ.reduction(b,BKZ.Param(block_size=28,max_loops=3),U=u)
 U=np.array([[int(u[i,j])for j in range(36)]for i in range(36)],dtype=np.int64);assert abs(U).max()<100000;assert np.array_equal(U@B,np.array([[int(b[i,j])for j in range(36)]for i in range(36)],dtype=np.int64))
 norm=np.einsum('vi,ij,vj->v',U,G,U);short=U[norm<8-1e-6];abc=float(np.prod(scale));center=(2**18/3**9)/(abc**6);row={'iteration':iteration,'scale':scale.tolist(),'abc':abc,'conditional_center':center,'basis_minimum':float(norm.min()),'cuts':len(short),'constraint_types':len(co)};rows.append(row);print(row,'seconds',time.time()-start,flush=True)
 if not len(short):
  np.savez(P/'CTturyn1_scale_candidate.npz',Q=G,U=U,scales=scale);np.savetxt(P/'CTturyn1_scale_candidate_R.txt',np.linalg.cholesky(U@G@U.T).T);break
 new=np.array([np.einsum('vi,ij,vj->v',short,p,short)for p in parts]).T;oldlen=len(co);co=np.unique(np.r_[co,new],axis=0)
 if len(co)==oldlen:raise RuntimeError('No progress from a rejecting vector')
(P/'CTturyn1_scale_cut.json').write_text(json.dumps(rows,indent=2)+'\n')
