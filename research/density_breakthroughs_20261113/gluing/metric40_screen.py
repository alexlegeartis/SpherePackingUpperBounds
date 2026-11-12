"""Density-specific metric screen of the P48/E8 40-dimensional lattice.
Uses exact min6 witnesses from affine-shell differences; numerical LP only.
"""
from pathlib import Path
import numpy as np,json,time,math,warnings,argparse
from scipy.optimize import linprog,OptimizeWarning
warnings.simplefilter('ignore',OptimizeWarning)
P=Path(__file__).resolve().parent;S=P.parents[3]/'KissingNumbers/verifications/improved/dim40-affine-shells/source'
ap=argparse.ArgumentParser();ap.add_argument('--contacts-only',action='store_true');args=ap.parse_args()
start=time.time();Q=np.array(json.loads((S/'gram.json').read_text()),dtype=np.int64);C=np.load(S/'core.npz')['coefficients'];CF=C.astype(np.float32);rng=np.random.default_rng(40);parts=[];anchors=rng.choice(len(C),20,replace=False)
for i,a in enumerate(anchors):
 z=C[a].astype(np.int64);w=Q@z;assert 40*int(np.max(abs(C)))*int(np.max(abs(w)))<2**24
 dots=CF@w.astype(np.float32);sel=np.flatnonzero(dots==27);diff=C[sel].astype(np.int16)-z;assert not np.any(diff%3);V=diff//3;assert np.all(np.einsum('bi,ij,bj->b',V,Q,V)==6);parts.append(V);print('anchor',i,'contacts',len(V),'seconds',time.time()-start,flush=True)
V=np.unique(np.concatenate(parts),axis=0);sign=np.sign(V[np.arange(len(V)),np.argmax(V!=0,axis=1)]);V=np.unique(V*sign[:,None],axis=0);np.savez(P/'contacts40.npz',Q=Q,V=V,anchors=anchors);print('unique',len(V),flush=True)
if args.contacts_only:raise SystemExit(0)
L=np.linalg.cholesky(Q);n=40;ii,jj=np.triu_indices(n);f=np.where(ii==jj,1.,2.);o=np.where(ii==jj,1.,0.);X=V@L;A=X[:,ii]*X[:,jj]*f
H=A.T@A;eig=np.linalg.eigvalsh(H);rank=np.sum(eig>eig[-1]*1e-10);print('rank',rank,'small_eig',eig[:8].tolist(),flush=True)
sel=rng.choice(len(V),min(3500,len(V)),replace=False);active=A[sel]
for trial in range(12):
 s=linprog(o,A_ub=-active,b_ub=np.zeros(len(active)),bounds=[(-1,1)]*len(ii),method='highs',options={'threads':1,'time_limit':55})
 print('LP',trial,s.message,None if s.x is None else s.fun,'time',time.time()-start,flush=True)
 if not s.success:break
 vals=A@s.x;bad=np.flatnonzero(vals<-1e-6)
 if not len(bad) or s.fun>-1e-7:break
 use=bad[np.argsort(vals[bad])[:300]];active=np.vstack([active,A[use]])
out={'dimension':40,'known_minimum':6,'known_Gram_determinant':6561,'center_density':3**16/2**20,'sampled_exact_contacts':len(V),'sampled_contact_rank':int(rank),'symmetric_space_dimension':len(ii),'LP_status':s.message,'LP_trace':float(s.fun) if s.fun is not None else None,'negative_sampled_constraints':len(bad) if s.success else None,'seconds':time.time()-start,'scope':'numerical metric screening over exact contact vectors, not global extremality or minimum certificate'}
if s.success:
 D=np.zeros((n,n));D[ii,jj]=s.x;D[jj,ii]=s.x;np.savez(P/'metric40_direction.npz',orthogonal_D=D,original_D=L@D@L.T)
(P/'metric40_results.json').write_text(json.dumps(out,indent=2));print(json.dumps(out),flush=True)
