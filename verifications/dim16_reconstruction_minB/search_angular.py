import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,numpy as np
from scipy.optimize import linprog,minimize
from pathlib import Path
from fractions import Fraction as F
from flint import fmpq,fmpq_poly
D=Path(__file__).resolve().parent
p=np.array(json.loads((D/'fit.json').read_text())['parameters']);VK=np.array([1.02,1.05,1.1,1.2,1.4,2]);LK=np.array([1.2,1.4,1.5,1.6,1.65,1.75]); DEG=19
def phi(x):
 return p[0]*(x-1)+p[7]*np.log(x)+np.maximum(x[...,None]-VK,0)@p[1:7]+np.maximum(np.log(x[...,None]/LK),0)@p[8:]
def G(x):
 x=np.asarray(x);a=[np.ones_like(x),x]
 for k in range(1,DEG):a.append(((2*k+14)*x*a[-1]-k*a[-2])/(k+14))
 return np.stack(a,axis=-1)
def rows(t,u):
 y=u*u;z=1+y-2*u*t;w=np.sqrt(z)
 return G(t)+np.exp(phi(y)-phi(z))[...,None]*G((1-u*t)/w)+np.exp(-phi(z))[...,None]*G((u-t)/w)
t,u=np.meshgrid(np.linspace(-1,.5,401),np.r_[np.linspace(1,1.5,151),np.linspace(1.5,4,101)])
m=2*t*u<=1+1e-14;t=t[m];u=u[m]
# Explicit sorted boundary
ub=np.linspace(1,4,1001);t=np.r_[t,1/(2*ub)];u=np.r_[u,ub]
x=np.linspace(-1,.5,10001);a=G(x)
for i in range(12):
 r=rows(t,u);A=np.vstack([r[:,1:],a[:,1:]]);b=np.r_[-r[:,0]-0.0002,-a[:,0]-0.0002]
 s=np.maximum(abs(A).max(axis=1),1e-10)
 fit=linprog(np.ones(DEG),A_ub=A/s[:,None],b_ub=b/s,bounds=(0,None),method='highs',options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
 assert fit.success,fit.message
 c=np.r_[1,fit.x];v=r@c;new=[]
 for j in np.argsort(v)[-50:]:
  opt=minimize(lambda q: -float(rows(np.array(q[0]),np.array(q[1]))@c),[t[j],u[j]],bounds=[(-1,.5),(1,4)],constraints=[{'type':'ineq','fun':lambda q:1-2*q[0]*q[1]}],method='SLSQP',options={'maxiter':100,'ftol':1e-11})
  if 2*opt.x[0]*opt.x[1]<=1+1e-10 and -opt.fun>-.00019:new.append(opt.x)
 print(i,1+fit.fun,max(v),len(new),flush=True)
 (D/'search_candidate.json').write_text(json.dumps({'c':c.tolist(),'B':sum(c),'certified':False},indent=2))
 if not new:break
 e=np.array(new);t=np.r_[t,e[:,0]];u=np.r_[u,e[:,1]]
cc=[fmpq(F(str(z)).numerator,F(str(z)).denominator) for z in c];X=fmpq_poly([0,1]);gg=[fmpq_poly([1]),X]
for k in range(1,DEG):gg.append(((2*k+14)*X*gg[-1]-k*gg[-2])/(k+14))
q=sum((a*b for a,b in zip(cc,gg)),fmpq_poly())
out={'dimension':16,'B_exact':str(sum(cc)),'gegenbauer_coefficients':list(map(str,cc)),'polynomial_coefficients':list(map(str,q.coeffs()))}
(D/'angular.json').write_text(json.dumps(out,indent=2)+'\n')
