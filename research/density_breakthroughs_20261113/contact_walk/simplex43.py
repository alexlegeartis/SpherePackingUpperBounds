"""Explore rank-five even Gram candidates for six-translate antipode packings.

Candidate Gram matrices are NOT embeddings and do not constitute packings.
"""
from pathlib import Path
import itertools,json,time
import numpy as np
P=Path(__file__).resolve().parent
n=5;i,j=np.triu_indices(n,1);base=226.5172748890234
z=np.array([a for a in itertools.product([-1,0,1],repeat=n) if any(a) and next(x for x in a if x)!=-1],dtype=float)
total=4**10;best=0.;records=[];tested=0;starts=time.time()
for start in range(0,total,16384):
 ids=np.arange(start,min(start+16384,total),dtype=np.int64)
 a=-((ids[:,None]//(4**np.arange(10)))%4).astype(float)
 K=np.broadcast_to(np.eye(n)*6,(len(ids),n,n)).copy();K[:,i,j]=a;K[:,j,i]=a
 ev=np.linalg.eigvalsh(K);ok=ev[:,0]>1e-7;K=K[ok];ids=ids[ok]
 if not len(K):continue
 H=np.linalg.inv(K);diag=np.diagonal(H,axis1=1,axis2=2)
 beta=np.maximum(diag.max(axis=1),(diag[:,:,None]+diag[:,None,:]-2*H).max(axis=(1,2)))
 det=np.linalg.det(K);positive=beta<6
 scores=np.zeros(len(K));scores[positive]=6*np.exp(21.5*np.log((6-beta[positive])/4)-.5*np.log(det[positive]))
 good=np.where(scores>max(best,base*.999))[0]
 if len(good):
  norms=np.einsum('zi,bij,zj->bz',z,K[good],z)
  for b,mi in zip(good,norms.min(axis=1)):
   if mi<6-1e-7:continue
   item={'index':int(ids[b]),'Gram':K[b].astype(int).tolist(),'det_float':float(det[b]),'beta_float':float(beta[b]),'center_density_float':float(scores[b]),'small_integer_min':float(mi),'embedded':False}
   if scores[b]>best+1e-8:
    best=float(scores[b]);records.append(item);print(json.dumps(item),flush=True)
 tested+=len(K)
 if start%(16384*8)==0:print('progress',start,total,'best',best,'seconds',time.time()-starts,flush=True)
(P/'simplex43.json').write_text(json.dumps({'box_offdiagonals':[-3,0],'diagonal':6,'total_candidates':total,'positive_definite':tested,'records':records,'scope':'finite abstract Gram and standard unimodular simplex search; no embedding or global bound'},indent=2)+'\n')
