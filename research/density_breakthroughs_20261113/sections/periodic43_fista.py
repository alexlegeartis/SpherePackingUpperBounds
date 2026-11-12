"""Independent periodic43 coordinate audit and a numerical positive-weight screen.

The finite exact contacts are valid independently of pool completeness.  A
numerical nonnegative solution to A.T w = trace is evidence against first-order
density descent; this program does not certify exact local optimality.
"""
from pathlib import Path
import json,re,time,sys
import numpy as np
import sympy as s
from scipy.linalg import null_space
from sympy.matrices.normalforms import smith_normal_form
P=Path(__file__).resolve().parent;BASE=P.parent
C=BASE/'seed_deformations/external_sphere_packing/packings_out/records/dim43-delta226/P48n-081e3830.cert.txt'
def matrix(text,key):
 a=text.index('[[',text.index(key));b=text.index(']]',a)+2
 return np.array([[int(x)for x in re.findall(r'-?\d+',row)]for row in text[a:b].splitlines()],dtype=np.int64)
start=time.time();text=C.read_text();G=matrix(text,'ambient Gram G');T=matrix(text,"simplex-basis coefficient rows T'")
K=T@G@T.T;assert s.Matrix(K).det()==3072
sm=smith_normal_form(s.Matrix(T),domain=s.ZZ);assert all(abs(sm[i,i])==1 for i in range(5))
adj=np.array(s.Matrix(K).inv()*3072,dtype=np.int64);assert np.array_equal(K@adj,3072*np.eye(5,dtype=np.int64))
pool=np.load(BASE/'gluing/P48n_glue_pool.npz');assert np.array_equal(pool['G'],G)
V=pool['V'];F=V.astype(np.int64)@G@T.T
S=np.vstack([np.zeros(5,dtype=np.int64),np.eye(5,dtype=np.int64)])
contacts=[];pairs=[];counts=[]
for i in range(6):
 for j in range(i+1,6):
  u=S[i]-S[j];b=int(u@adj@u);assert b in (768,896)
  if b!=896:continue
  a=V[np.all(F==u,axis=1)];b=-V[np.all(F==-u,axis=1)]
  W=np.concatenate([a,b]).astype(np.int64)
  assert np.all(np.einsum('ij,ij->i',W@G,W)==6)
  contacts.extend(W);pairs.extend([(i,j)]*len(W));counts.append([i,j,len(W)])
V=np.array(contacts,dtype=np.int64);pairs=np.array(pairs,dtype=np.int64)
full_count=len(V);rng=np.random.default_rng(43911)
count=min(full_count,int(sys.argv[1])if len(sys.argv)>1 else 12000)
ids=rng.choice(full_count,count,replace=False);V=V[ids];pairs=pairs[ids]
L=np.linalg.cholesky(G);E=null_space(T@L);W=V@L@E
assert np.max(abs(np.einsum('ij,ij->i',W,W)-137/24))<1e-8
assert E.shape==(48,43) and np.linalg.norm((T@L)@E)<1e-10
n=43;ii,jj=np.triu_indices(n);q=len(ii)
A=np.zeros((count,q+5*n));A[:,:q]=W[:,ii]*W[:,jj]*(1+(ii!=jj))
for i in range(1,6):A[:,q+(i-1)*n:q+i*n]=2*W*((pairs[:,0]==i).astype(int)-(pairs[:,1]==i).astype(int))[:,None]
goal=np.zeros(A.shape[1]);goal[:q]=(ii==jj)
# Audit the symmetric-metric convention and oriented coset translation signs.
test=rng.normal(size=A.shape[1]);D=np.zeros((n,n));D[ii,jj]=test[:q];D[jj,ii]=test[:q]
shifts=np.vstack([np.zeros(n),test[q:].reshape(5,n)])
direct=np.einsum('ij,ij->i',W@D,W)+2*np.einsum('ij,ij->i',W,shifts[pairs[:,0]]-shifts[pairs[:,1]])
assert np.max(abs(direct-A@test))<1e-9
eig,U=np.linalg.eigh(A.T@A);keep=eig>eig[-1]*1e-10
t=U[:,keep]@((U[:,keep].T@goal)/eig[keep]);weights=A@t
report={'contacts_available':full_count,'contacts_used':count,'pair_counts':counts,'metric_variables':q,'translate_variables':5*n,'total_variables':A.shape[1],'numerical_rank':int(keep.sum()),'smallest_eigenvalues':eig[:12].tolist(),'least_norm_min_weight':float(weights.min()),'least_norm_negative_weights':int((weights<0).sum()),'least_norm_residual':float(np.linalg.norm(A.T@weights-goal)),'normal_gram_determinant':3072,'normal_basis_saturated':True,'contact_squared_distance':'137/24','seed_center_density':float(6*(137/96)**21.5/np.sqrt(3072)),'scope':'exact individual ambient contact norms and saturation; numerical full periodic first-order screen on a finite subset; no new packing or exact local-optimality claim'}
print(json.dumps(report),flush=True)
if weights.min()<-1e-10:
 X=(A@U[:,keep])/np.sqrt(eig[keep]);target=(U[:,keep].T@goal)/np.sqrt(eig[keep])
 best=np.maximum(weights,0);y=best.copy();t=1.;began=time.time()
 for it in range(20001):
  residual=X.T@y-target;nxt=np.maximum(0,y-X@residual)
  nt=(1+np.sqrt(1+4*t*t))/2;y=nxt+(t-1)/nt*(nxt-best);best=nxt;t=nt
  if it%100==0:
   err=float(np.linalg.norm(X.T@best-target))
   if it%1000==0:print('FISTA',it,err,'seconds',time.time()-began,flush=True)
   if err<1e-12 or time.time()-began>120:break
 weights=best
 report.update(fista_iterations=it,fista_transformed_residual=err)
report.update(final_min_weight=float(weights.min()),final_positive_weights=int((weights>1e-12).sum()),final_goal_residual=float(np.linalg.norm(A.T@weights-goal)),elapsed=time.time()-start)
assert np.max(abs(V))<=32767
np.savez_compressed(P/'periodic43_fista.npz',G=G,T=T,V=V.astype(np.int16),pairs=pairs.astype(np.uint8),E=E,weights=weights)
(P/'periodic43_fista_report.json').write_text(json.dumps(report,indent=2));print('FINAL',json.dumps(report),flush=True)
