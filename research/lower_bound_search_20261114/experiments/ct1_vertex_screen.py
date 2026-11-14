"""Floating discovery in ten exact contact-preserving CT vertex directions."""
from pathlib import Path
import json,time
import numpy as np
from scipy.optimize import linprog
from fpylll import IntegerMatrix,BKZ
P=Path(__file__).resolve().parent;start=time.time();x=np.load(P/'CTturyn1_vertex_flex.npz');Q=x['Q'];D=x['D'];den=int(x['denom']);V=np.loadtxt(P/'CTturyn1_vertex_short10.txt',dtype=np.int64);first=np.argmax(V!=0,axis=1);V=V[V[np.arange(len(V)),first]>0];norm=np.einsum('vi,ij,vj->v',V,Q,V);assert norm.min()==128
i,j=np.triu_indices(36);K=(D[:,i,j]*np.where(i==j,1,2)).T.astype(float);rows=[]
for begin in range(0,len(V),10000):
 W=V[begin:begin+10000];A=np.rint((W[:,i]*W[:,j]).astype(float)@K).astype(np.int64);rows.append(np.c_[norm[begin:begin+len(W)],A])
constraints=np.unique(np.concatenate(rows),axis=0);constraints=constraints[np.any(constraints[:,1:]!=0,axis=1)];print('constraints',len(constraints),'seconds',time.time()-start,flush=True);rng=np.random.default_rng(20261114);history=[]
for cutting_round in range(15):
 candidates={}
 for trial in range(180):
  c=rng.normal(size=len(D))
  for inner in range(3):
   result=linprog(c,A_ub=-constraints[:,1:],b_ub=constraints[:,0]-128,bounds=[(-20,20)]*len(D),method='highs')
   if not result.success:break
   t=result.x;G=Q+np.einsum('a,aij->ij',t,D);evals=np.linalg.eigvalsh(G)
   if evals[0]<1e-7:break
   delta=2**18/np.exp(np.linalg.slogdet(G/den)[1]/2);key=tuple(np.round(t,9));candidates[key]=(delta,t,G)
   inverse=np.linalg.inv(G);c=np.einsum('ij,aji->a',inverse,D)
 print('round',cutting_round,'positive candidates',len(candidates),'seconds',time.time()-start,flush=True)
 if not candidates:break
 ordered=sorted(candidates.values(),key=lambda v:-v[0]);newcuts=[]
 for delta,t,G in ordered[:12]:
  B=np.rint(np.linalg.cholesky(G/den)*1e7).astype(np.int64);b=IntegerMatrix.from_matrix(B.tolist());u=IntegerMatrix.identity(36);BKZ.reduction(b,BKZ.Param(block_size=28,max_loops=3),U=u);U=np.array([[int(u[i,j])for j in range(36)]for i in range(36)],dtype=np.int64);assert np.array_equal(U@B,np.array([[int(b[i,j])for j in range(36)]for i in range(36)],dtype=np.int64));norms=np.einsum('vi,ij,vj->v',U,G,U);W=U[norms<128-1e-6];row={'round':cutting_round,'conditional_center_density':float(delta),'t':t.tolist(),'basis_minimum':float(norms.min()/den),'rejecting_rows':len(W)};history.append(row);print(row,'seconds',time.time()-start,flush=True)
  if not len(W):
   np.savez(P/f'CTturyn1_flex_candidate_{cutting_round}.npz',Q=G/den,U=U,t=t,D=D,base=Q,denom=den);np.savetxt(P/f'CTturyn1_flex_candidate_{cutting_round}_R.txt',np.linalg.cholesky(U@G@U.T/den).T)
  else:
   A=np.array([np.einsum('vi,ij,vj->v',W,d,W)for d in D]).T;N=np.einsum('vi,ij,vj->v',W,Q,W);newcuts.extend(np.c_[N,A])
 (P/'CTturyn1_vertex_screen.json').write_text(json.dumps(history,indent=2)+'\n')
 if not newcuts:break
 constraints=np.unique(np.r_[constraints,newcuts],axis=0)
 np.savez(P/'CTturyn1_vertex_constraints.npz',constraints=constraints,Q=Q,D=D,denom=den)
