"""Exact contact-preserving directions for the alternative CT scalar vertex."""
from pathlib import Path
import numpy as np,time
from scipy.linalg import qr
from flint import fmpz_mat
P=Path(__file__).resolve().parent;start=time.time();x=np.load(P/'CTturyn1_scale_vertex.npz');Q=np.rint(x['Q']*16).astype(np.int64);V=np.loadtxt(P/'CTturyn1_scale_vertex_vectors.txt',dtype=np.int64);V=np.unique(np.array([v if next((int(z)for z in v if z),0)>0 else-v for v in V]),axis=0);assert np.all(np.einsum('vi,ij,vj->v',V,Q,V)==128)
i,j=np.triu_indices(36);A=V[:,i]*V[:,j]*np.where(i==j,1,2);_,R,p=qr(A.T.astype(float),mode='economic',pivoting=True);rank=int(np.sum(abs(np.diag(R))>1e-7));print('rank',rank,'QR seconds',time.time()-start,flush=True);M=fmpz_mat(A[p[:rank]].tolist());K,nullity=M.nullspace();assert nullity==666-rank;print('exact nullity',nullity,'seconds',time.time()-start,flush=True);L=fmpz_mat([[int(K[i,j])for i in range(666)]for j in range(nullity)]).lll();B=np.array(L.tolist(),dtype=np.int64);assert not (A@B.T).any();D=np.zeros((nullity,36,36),dtype=np.int64);D[:,i,j]=B;D[:,j,i]=B;print('D max',np.max(abs(D)),'seconds',time.time()-start,flush=True);np.savez(P/'CTturyn1_vertex_flex.npz',Q=Q,D=D,V=V,denom=16)
