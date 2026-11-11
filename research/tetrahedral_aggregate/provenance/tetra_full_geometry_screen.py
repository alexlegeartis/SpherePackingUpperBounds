"""Bounded numerical tetrahedron LP. No continuum sign/density certificate."""
import itertools,json,time
from pathlib import Path
import numpy as np
from scipy.special import eval_gegenbauer,poch,comb,factorial
from scipy.optimize import linprog,minimize
ROOT=Path(__file__).resolve().parent

def far_bernstein(d,D,hi):
 lam=(d-2)/2;parts=[]
 for left,right in zip(np.linspace(-1,hi,41)[:-1],np.linspace(-1,hi,41)[1:]):
  jets=np.zeros((D+1,D+1))
  for j in range(D+1):
   for ell in range(j,D+1):jets[j,ell]=2**j*poch(lam,j)*eval_gegenbauer(ell-j,lam+j,left)/eval_gegenbauer(ell,lam,1)/factorial(j)
  B=np.zeros((D+1,D+1))
  for i in range(D+1):
   for j in range(i+1):B[i,j]=comb(i,j)/comb(D,j)*(right-left)**j
  parts.append(B@jets)
 M=np.vstack(parts);return M/np.maximum(np.max(np.abs(M),axis=1,keepdims=True),1e-200)

EDGES=list(itertools.combinations(range(4),2));INDEX={x:i for i,x in enumerate(EDGES)}
def edge(X,i,j):return X[...,INDEX[tuple(sorted((i,j)))]]
def cosines(X):
 X=np.asarray(X);answer=[]
 for i,j in EDGES:
  k,l=[z for z in range(4)if z not in(i,j)]
  a=edge(X,i,j);u=edge(X,i,k);v=edge(X,j,k);w=edge(X,i,l);z=edge(X,j,l);c=edge(X,k,l)
  R=(u+v)/2-a/4-(u-v)**2/(4*a);S=(w+z)/2-a/4-(w-z)**2/(4*a)
  delta=(u-v-w+z)**2/(4*a)
  answer.append((R+S+delta-c)/(2*np.sqrt(R*S)))
 return np.stack(answer,axis=-1)
def rows(X,n,D):
 t=cosines(X);lam=(n-3)/2
 return np.stack([np.sum(eval_gegenbauer(k,lam,t)/eval_gegenbauer(k,lam,1),axis=-1)for k in range(D+1)],axis=-1)/6

result=[];rng=np.random.default_rng(20261111);start=time.time()
for n in [6,16,22]:
 for q in [1.05,1.1,1.2]:
  D=24;far=1-q/(2*(q-.25));A=far_bernstein(n-1,D,far)
  X=np.r_[np.asarray(list(itertools.product(np.linspace(1,q,3),repeat=6))),rng.uniform(1,q,(2500,6))]
  M=rows(X,n,D);history=[]
  for rnd in range(3):
   MM=np.r_[A,M]
   lp=linprog(np.ones(D),A_ub=MM[:,1:],b_ub=-MM[:,0]-1e-8,bounds=[(0,None)]*D,method='highs')
   if not lp.success:break
   cs=np.r_[1.,lp.x];peaks=[];values=M@cs
   for x in X[np.argsort(values)[-12:]]:
    opt=minimize(lambda z:-float(rows(z,n,D)@cs),x,bounds=[(1,q)]*6,method='L-BFGS-B',options={'ftol':1e-13,'gtol':1e-10,'maxiter':80})
    peaks.append((float(-opt.fun),opt.x))
   extra=np.asarray([x for v,x in peaks if v>-1e-9]);worst=max(v for v,x in peaks)
   info={'round':rnd,'B_numerical':float(sum(cs)),'worst_local_tetra':worst,'added':len(extra)};history.append(info)
   if len(extra):X=np.r_[X,extra];M=np.r_[M,rows(extra,n,D)]
   else:break
  rec={'dimension':n,'q':q,'history':history,'coefficients':cs.tolist(),'certified':False,'warning':'Sampled six-dimensional tetrahedron signs; no global proof.'}
  result.append(rec);print(n,q,history,flush=True)
  (ROOT/'tetra_full_geometry_screen.json').write_text(json.dumps(result,indent=2)+'\n')
print('seconds',time.time()-start,flush=True)
