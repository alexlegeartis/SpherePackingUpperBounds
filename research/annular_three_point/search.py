"""Numerical discovery only: annular Cohn-de Laat-Salmon three-point SDP.
No result from this script is a proof. Export is intentionally marked uncertified.
"""
import argparse,json,math,time
from pathlib import Path
import numpy as np
from scipy.special import eval_genlaguerre,gammaln,eval_gegenbauer
import cvxpy as cp

def solve(n,degree=39,nodes=5,angular=16,ratio=1.2,solver='CLARABEL'):
 if not 3<=n<=48: raise ValueError("Supported dimensions are 3 through 48")
 if degree<1 or nodes<2 or angular<1 or not 1<ratio<2: raise ValueError("Invalid search parameters")
 a=n/2; r2=n/8; R2=r2*ratio; u=np.linspace(r2,R2,nodes); scale=np.exp(-np.pi*u/2)
 coeff=cp.Variable(degree+1); blocks=[cp.Variable((nodes+1,nodes+1),symmetric=True)]+[cp.Variable((nodes,nodes),symmetric=True) for _ in range(angular)]
 norms=np.exp(gammaln(np.arange(degree+1)+a)-gammaln(a)-gammaln(np.arange(degree+1)+1))
 def lag(x):return np.array([eval_genlaguerre(k,a-1,2*np.asarray(x))/norms[k] for k in range(degree+1)]).T
 def weights(v):
  z=(v-r2)/(R2-r2)*(nodes-1); i=min(nodes-2,max(0,int(z))); w=np.zeros(nodes); w[i]=1-z+i; w[i+1]=z-i; return w
 cons=[B>>1e-5*np.eye(B.shape[0]) for B in blocks]
 cons += [np.power(-1.,np.arange(degree+1))@coeff==1]
 # Sampling only. Exact polynomial tail and full angular/annular certification
 # are separate gates and are not implied by these constraints.
 xs=np.unique(np.r_[np.linspace(0,80,500),np.geomspace(80,5000,100)])
 for sign,xx in [(1,xs),(-1,np.unique(np.r_[np.linspace(np.pi*R2,80,500),np.geomspace(80,5000,100)]))]:
  mat=lag(xx); mat/=np.maximum(np.linalg.norm(mat,axis=1,keepdims=True),1)
  if sign==1:cons += [mat@cp.multiply(np.power(-1.,np.arange(degree+1)),coeff)>=1e-7]
  else:cons += [mat@coeff<=-1e-7]
 for v in np.linspace(r2,R2,8*(nodes-1)+1):
  w=weights(v)*scale*np.exp(np.pi*v/2)
  diag=sum(cp.sum(cp.multiply(np.outer(w,w),B[1:,1:] if j==0 else B)) for j,B in enumerate(blocks))
  cross=2*np.exp(np.pi*v/2)*(w@blocks[0][1:,0])
  cons += [lag(np.pi*v)@coeff+cross+diag<=-1e-5]
 fine=np.linspace(r2,R2,2*(nodes-1)+1)
 for i,v in enumerate(fine):
  for z in fine[i:]:
   w=weights(v)*scale; ww=weights(z)*scale; outer=np.outer(w,ww)/(np.exp(-np.pi*(v+z)/2))
   cs=cp.hstack([cp.sum(cp.multiply(outer,B[1:,1:] if j==0 else B)) for j,B in enumerate(blocks)])
   t=np.linspace(-1,(v+z-r2)/(2*np.sqrt(v*z)),60)
   g=np.array([eval_gegenbauer(j,(n-2)/2,t)/eval_gegenbauer(j,(n-2)/2,1.) for j in range(angular+1)]).T
   cons += [g@cs<=-1e-5]
 objective=cp.sum(coeff)+blocks[0][0,0]
 p=cp.Problem(cp.Minimize(objective),cons)
 opts={'max_iter':300, 'tol_gap_abs':1e-10, 'tol_feas':1e-10, 'tol_gap_rel':1e-10} if solver=='CLARABEL' else {'eps':1e-6,'max_iters':10000}
 start=time.time();p.solve(solver=solver,**opts)
 result={'dimension':n,'status':str(p.status),'certified':False,'warning':'Sampled SDP output; not a rigorous bound','seconds':time.time()-start,'radial_degree':degree,'radial_nodes':nodes,'angular_degree':angular,'exclusion_squared':str(r2),'truncation_squared':str(R2)}
 if coeff.value is not None:
  sc=np.r_[1,scale]; mats=[blocks[0].value*np.outer(sc,sc)]+[B.value*np.outer(scale,scale) for B in blocks[1:]]
  result.update(coefficients=coeff.value.tolist(),matrices=[B.tolist() for B in mats],objective=float(p.value),density_candidate=float(p.value)*np.pi**a/math.gamma(a+1)*(r2/4)**a,max_constraint_violation=max(float(np.max(c.violation())) for c in cons))
 return result
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--dimensions',type=int,nargs='+',default=[17,18,22,23]);ap.add_argument('--degree',type=int,default=39);ap.add_argument('--nodes',type=int,default=5);ap.add_argument('--angular',type=int,default=16);ap.add_argument('--solver',default='CLARABEL');ap.add_argument('--output',default='candidates.json');args=ap.parse_args();out=[]
 for n in args.dimensions:
  try:r=solve(n,args.degree,args.nodes,args.angular,solver=args.solver)
  except Exception as e:r={'dimension':n,'certified':False,'error':str(e)}
  out.append(r);Path(args.output).write_text(json.dumps(out,indent=2));print({k:v for k,v in r.items() if k not in ['coefficients','matrices']},flush=True)
