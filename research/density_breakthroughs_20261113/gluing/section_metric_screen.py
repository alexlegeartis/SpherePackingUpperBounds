"""Numerical contact-cone screens for the reference 39D and 40D lattices."""
from pathlib import Path
import numpy as np,sympy as sp,json,time,warnings,argparse
from scipy.optimize import linprog,OptimizeWarning
warnings.simplefilter('ignore',OptimizeWarning)
P=Path(__file__).resolve().parent;ap=argparse.ArgumentParser();ap.add_argument('--dimension',type=int,choices=(39,40),default=39);ap.add_argument('--seconds',type=float,default=100);a=ap.parse_args();n=a.dimension
j=np.load(P/'contacts40.npz');Q=j['Q'];V=j['V'];J=np.array(sp.Matrix(Q).inv()*3,dtype=np.int64);eye=np.eye(40,dtype=np.int64);short=[]
for i in range(40):
 if J[i,i]==10:short.append(eye[i])
 for k in range(i):
  for s in(-1,1):
   if J[i,i]+J[k,k]+2*s*J[i,k]==10:short.append(eye[i]+s*eye[k])
if n==40:C=np.zeros((0,40),dtype=np.int64);K=eye
elif n==39:C=np.array([short[2]])
if n<40:
 R,piv=sp.Matrix(C).rref();free=[i for i in range(40)if i not in piv];K=np.zeros((n,40),dtype=np.int64);K[:,free]=np.eye(n,dtype=np.int64)
 for i,p in enumerate(piv):K[:,p]=-np.array(R[i,free],dtype=np.int64).reshape(-1)
 assert not np.any(K@C.T);keep=np.flatnonzero(np.all(V@C.T==0,axis=1));V=V[keep][:,free]
G=K@Q@K.T;det=sp.Matrix(G).det();assert np.all(np.einsum('bi,ij,bj->b',V,G,V)==6);den=float((sp.Rational(3,2))**sp.Rational(n,2)/sp.sqrt(det));print('seed',n,'det',det,'density',den,'contacts',len(V),flush=True)
L=np.linalg.cholesky(G);ii,jj=np.triu_indices(n);factor=np.where(ii==jj,1.,2.);o=np.where(ii==jj,1.,0.);coords=V@L;A=coords[:,ii]*coords[:,jj]*factor
st=time.time();e=np.linalg.eigvalsh(A.T@A);rank=int(np.sum(e>e[-1]*1e-10));print('rank',rank,len(ii),'eigen',e[:5].tolist(),flush=True)
rng=np.random.default_rng(n);ids=rng.choice(len(V),min(3000,len(V)),replace=False);act=A[ids]
for k in range(10):
 res=linprog(o,A_ub=-act,b_ub=np.zeros(len(act)),bounds=[(-1,1)]*len(o),method='highs-ipm',options={'threads':1,'time_limit':a.seconds,'ipm_optimality_tolerance':1e-9})
 print('LP',k,res.message,'fun',res.fun,'time',time.time()-st,flush=True)
 if not res.success:break
 vals=A@res.x;bad=np.flatnonzero(vals<-1e-6)
 if not len(bad)or res.fun>-1e-6:break
 ids=bad[np.argsort(vals[bad])[:400]];act=np.vstack([act,A[ids]])
out={'dimension':n,'seed_Gram_determinant':str(det),'seed_minimum':6,'seed_center_density':den,'exact_contact_count':len(V),'numerical_contact_rank':rank,'symmetric_dimension':len(ii),'LP_status':res.message,'LP_trace':None if res.fun is None else float(res.fun),'negative_sampled_constraints':None if not res.success else int(len(bad)),'seconds':time.time()-st,'scope':'numerical contact cone screening; no record or global minimum claim'}
(P/f'section_metric{n}.json').write_text(json.dumps(out,indent=2));np.savez(P/f'section_metric{n}.npz',G=G,K=K,V=V,normals=C)
print(json.dumps(out),flush=True)
