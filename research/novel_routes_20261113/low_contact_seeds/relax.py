"""Experimental full-metric search from sparse-contact odd unimodular seeds.

Every minimum calculation here is a floating discovery calculation. A candidate
would need rational reconstruction and an independent minimum certificate.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import json,subprocess,time,sys
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from fpylll import IntegerMatrix,LLL
P=Path(__file__).resolve().parent
name=sys.argv[1];iterations=int(sys.argv[2]) if len(sys.argv)>2 else 60
out=P/name;out.mkdir(exist_ok=True)
Q=np.load(P/(name+'.npy')).astype(float);n=len(Q)
ii,jj=np.triu_indices(n);mult=np.where(ii==jj,1.,2.)
objective=(ii==jj).astype(float)
history=[];t0=time.time()
def shell(G,bound):
 np.savetxt(out/'Q.txt',G,fmt='%.17g');np.savetxt(out/'R.txt',np.linalg.cholesky(G).T,fmt='%.17g')
 cp=subprocess.run([str(P/'float_shell'),str(n),str(out/'Q.txt'),str(out/'R.txt'),str(bound),'150000000',str(out/'shell.txt'),'100000'],text=True,capture_output=True,check=True)
 stat=json.loads(cp.stdout)
 if stat['leaves']:
  v=np.loadtxt(out/'shell.txt',ndmin=2)[:,1:].astype(np.int64)
 else:v=np.zeros((0,n),np.int64)
 return v,stat
for iteration in range(iterations):
 B=np.linalg.cholesky(Q)
 Z=IntegerMatrix.from_matrix(np.rint(B*1e8).astype(np.int64).tolist());U=IntegerMatrix.identity(n)
 LLL.reduction(Z,U,delta=.99)
 T=np.array([[int(U[a,b])for b in range(n)]for a in range(n)],np.int64)
 Q=T@Q@T.T;Q=(Q+Q.T)/2
 B=np.linalg.cholesky(Q);V,initial=shell(Q,3.15)
 if not len(V):raise RuntimeError('no shell vectors')
 minold=float(np.einsum('bi,ij,bj->b',V,Q,V).min())
 if minold<3-1e-6:
  Q*=3/minold;B=np.linalg.cholesky(Q)
 oldlog=float(np.linalg.slogdet(Q)[1]);accepted=False
 for cutround in range(18):
  W=V@B;norms=np.sum(W*W,axis=1);A=W[:,ii]*W[:,jj]*mult
  lp=linprog(objective,A_ub=-A,b_ub=norms-3,bounds=(-.25/n,.25/n),method='highs',options={'time_limit':25,'threads':1})
  if lp.x is None:break
  E=np.zeros((n,n));E[ii,jj]=lp.x;E[jj,ii]=lp.x
  proposal=B@(np.eye(n)+E)@B.T;proposal=(proposal+proposal.T)/2
  extra,stat=shell(proposal,3.000001)
  if len(extra):
   vals=np.einsum('bi,ij,bj->b',extra,proposal,extra);bad=extra[vals<3-2e-7]
  else:bad=np.zeros((0,n),np.int64)
  if len(bad):
   V=np.unique(np.vstack([V,bad]),axis=0);continue
  Q=proposal;accepted=True;break
 logdet=float(np.linalg.slogdet(Q)[1]);delta=float((3/4)**(n/2)*np.exp(-logdet/2))
 record={'iteration':iteration,'initial_shell':initial,'constraint_vectors':len(V),'cut_rounds':cutround+1,'accepted':accepted,'trace_objective':None if lp.x is None else float(lp.fun),'logdet':logdet,'conditional_density':delta,'seconds':time.time()-t0,'status':'floating discovery, not a certified packing'}
 history.append(record);print(json.dumps(record),flush=True)
 np.save(out/'current_gram.npy',Q);(out/'history.json').write_text(json.dumps(history,indent=2)+'\n')
 if not accepted or oldlog-logdet<1e-9:break
