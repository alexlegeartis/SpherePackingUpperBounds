"""Six-parameter contact-preserving Gram discovery, with exact witness cuts.
Floating oracle incompleteness never certifies a minimum.
"""
import os,json,pathlib,time,subprocess,sys
from fractions import Fraction
from math import lcm
import numpy as np
from scipy.optimize import linprog
from fpylll import IntegerMatrix,BKZ,LLL
p=pathlib.Path(__file__).resolve().parent/'dim48mod3';d=np.load(p/'flexes.npz');H=d['H'];Ds=d['Ds'];V=d['V'];n=len(H);D=Ds/206.;base_detlog=np.linalg.slogdet(H)[1];record=(1.5)**24
st=np.load(p/'short_reduced.npz');T=st['T'];rng=np.random.default_rng(44)
# all cuts are integer v*D_i*v >= 206(8-v*H*v)
cache=p/'neighbor_cuts.npz'
if cache.exists():
 dat=np.load(cache);raw=dat['raw']
else:
 W=V[:6000];WV=W@H@T.T;qm=np.diag(T@H@T.T)
 aa=np.array([2*W@z@T.T for z in Ds]);dd=np.array([np.diag(T@z@T.T)for z in Ds]);qq=8+qm[None,:]+2*WV
 out=[]
 for sign in [1,-1]:
  qq=8+qm[None,:]+sign*2*WV;aa2=sign*aa+dd[:,None,:]
  keep=qq<=20
  out.append(np.column_stack([206*(8-qq[keep]),aa2[:,keep].T]))
 raw=np.unique(np.concatenate(out).astype(np.int64),axis=0);raw=raw[np.any(raw[:,1:],axis=1)];np.savez_compressed(cache,raw=raw)
print('initialcuts',len(raw),flush=True)
extra=[];history=[];best=0.;start=time.time()

def cuts(vs):
 q=np.einsum('ij,jk,ik->i',vs,H,vs);a=np.array([np.einsum('ij,jk,ik->i',vs,z,vs)for z in Ds]).T
 return np.column_stack([206*(8-q),a]).astype(np.int64)

def metric(x):return H+np.einsum('i,ijk->jk',x,D)

def score(x):
 G=metric(x);e=np.linalg.eigvalsh(G)
 return np.exp(24*np.log(2)-np.linalg.slogdet(G)[1]/2) if e[0]>0 else 0.

def lpwalk(x,bound=3.0):
 global raw
 for it in range(30):
  Q=metric(x);gi=np.linalg.inv(Q);grad=np.einsum('ij,kji->k',gi,D)
  sc=np.maximum(1,np.max(abs(raw),axis=1));res=linprog(grad,A_ub=-raw[:,1:]/sc[:,None],b_ub=-raw[:,0]/sc,bounds=[(-bound,bound)]*6,method='highs',options={'threads':1})
  if not res.success: print('LPfailure',res.message,flush=True);return x
  y=res.x
  Qy=metric(y);e,E=np.linalg.eigh(Qy)
  if e[0]<=1e-7:
   vec=np.rint(E[:,0]*10).astype(np.int64)
   if not np.any(vec):vec=np.rint(E[:,0]*100).astype(np.int64)
   raw=np.unique(np.vstack([raw,cuts(vec[None,:])]),axis=0);continue
  if np.linalg.norm(y-x)<1e-7:return y
  x=y
 return x

def oracle(x,trial,cap=50000000):
 rational=[Fraction(float(t)).limit_denominator(1000000) for t in x]
 common=lcm(*(t.denominator for t in rational))
 if common>1000000000:raise ValueError('Rational candidate denominator too large for this int64 discovery implementation')
 coeff=np.array([int(t*common) for t in rational],dtype=np.int64);den=206*common
 if np.any(raw[:,1:].astype(object)@coeff.astype(object)<raw[:,0].astype(object)*common):
  raise ValueError('Rational candidate violates a stored exact cut; refine reconstruction before enumerating')
 Q=den*H+np.einsum('i,ijk->jk',coeff,Ds);Q=(Q+Q.T)//2
 W=np.vstack([np.eye(n,dtype=np.int64),V[rng.choice(len(V),200,replace=False)]])
 B0=np.rint(np.linalg.cholesky(Q/den)*1e7).astype(np.int64);BB=IntegerMatrix.from_matrix((W@B0).tolist());U=IntegerMatrix.identity(len(W));LLL.reduction(BB,delta=.99,U=U)
 U=np.array([[int(U[i,j])for j in range(len(W))]for i in range(len(W))],dtype=np.int64);R=(U@W)[-n:]
 BB=IntegerMatrix.from_matrix((R@B0).tolist());U=IntegerMatrix.identity(n);BKZ.reduction(BB,BKZ.Param(block_size=28,max_loops=2),U=U,float_type='long double');U=np.array([[int(U[i,j])for j in range(n)]for i in range(n)],dtype=np.int64);R=U@R;Qr=R@Q@R.T
 np.savetxt(p/'Qcand.txt',Qr,fmt='%d');np.savetxt(p/'Rcand.txt',np.linalg.cholesky(Qr).T,fmt='%.17g')
 run=subprocess.run([str(pathlib.Path(__file__).resolve().parent.parent/'bin'/'shell_discovery'),str(n),str(p/'Qcand.txt'),str(p/'Rcand.txt'),str(8*den-1),str(cap),str(p/'bad.txt'),'0'],capture_output=True,text=True)
 rec=json.loads(run.stdout);rec.update({'trial':trial,'coefficients_numerators':coeff.tolist(),'coefficients_denominator':common,'density_if_minimum8':score(coeff/common),'exact_candidate':True,'verification':'discovery only'})
 bad=p/'bad.txt'
 if bad.stat().st_size:
  vs=np.loadtxt(bad,dtype=np.int64,ndmin=2)[:,1:]@R
  q=np.einsum('ij,jk,ik->i',vs,Q,vs);assert np.all(q<8*den);rec['found_exact_minimum_upper']=float(q.min()/den);return vs,rec
 np.savez(p/f'candidate_{trial}.npz',Q=Q,denominator=den,coefficients_numerators=coeff,coefficients_denominator=common,R=R,Qr=Qr)
 return np.empty((0,n),dtype=np.int64),rec

x=np.zeros(6)
for trial in range(50):
 x=lpwalk(x)
 s=score(x);print('candidate',trial,'density',s,'x',x.tolist(),'cuts',len(raw),'seconds',time.time()-start,flush=True)
 if not s:break
 vs,rec=oracle(x,trial)
 print(json.dumps(rec),flush=True);history.append(rec);(p/'descent_history.json').write_text(json.dumps(history,indent=2))
 if len(vs):
  raw=np.unique(np.vstack([raw,cuts(vs)]),axis=0);np.savez_compressed(p/'descent_cuts.npz',raw=raw)
  x=x*.98
 else:
  if s>best:
   best=s;np.savez(p/'best_candidate.npz',x=x,H=H,Ds=Ds,minimum_status='not certified')
  # Broad additional oracle basis to detect missed cuts.
  vs,rec=oracle(x,trial+1000,200000000)
  print(json.dumps(rec),flush=True);history.append(rec);(p/'descent_history.json').write_text(json.dumps(history,indent=2))
  if len(vs):raw=np.unique(np.vstack([raw,cuts(vs)]),axis=0);x=x*.98
  else:break
