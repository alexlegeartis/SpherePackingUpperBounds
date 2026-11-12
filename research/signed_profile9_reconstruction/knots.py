import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import experiment as base
import log_hinge as log
from basis import basis,old_to_new,new_to_old,knot_jacobian
D=Path(__file__).parent;j=json.loads((D/'fit9.json').read_text());K,q=old_to_new(j['parameters'],j['v_knots'],j['log_knots']);assert min(q)>-1e-10;q=np.maximum(q,0);K0=K.copy();rad=.4*np.minimum(K-np.r_[1,K[:-1]],np.r_[K[1:],2.01]-K);rad[-1]=0
base.KNOTS=K;log.ALL_KNOTS=K;log.LOG_KNOTS=K
class E(log.LogHingeExperiment):
 @staticmethod
 def basis(v):return basis(v,K)
e=E(9);raw=json.loads((D/'refit9.json').read_text());du=np.load(D/'duals9.npz');x=du['x'];v=x/(np.pi*9/8);h=np.exp(-basis(v,K)@q);mult=-du['dual']*raw['C']*np.exp(x)
grad=(mult*2*h*(1-h)/e.B)@knot_jacobian(v,q,K)*rad;grad/=np.linalg.norm(grad)
allr=e.grid(dense=True);vals=e.expression(allr,q);rng=np.random.default_rng(91);ids=np.unique(np.r_[np.argsort(vals)[-8000:],rng.choice(len(vals),2000,replace=False)]);x0=np.zeros(len(K))
for it in range(6):
 r={k:v[ids] for k,v in allr.items()};scale=np.maximum(1e-3,abs(r['q0'])+abs(r['q1'])+abs(r['q2']))
 def con(x,jac=False):
  kk=K0+x*rad;py=basis(r['y'],kk)@q;pz=basis(r['z'],kk)@q;v1=np.exp(py-pz)*r['q1'];v2=np.exp(-pz)*r['q2'];val=r['q0']+v1+v2
  if jac:
   jy=knot_jacobian(r['y'],q,kk);jz=knot_jacobian(r['z'],q,kk);return -(v1[:,None]*(jy-jz)-v2[:,None]*jz)*rad/scale[:,None]
  return (-val-1e-5)/scale
 fit=minimize(lambda x:float(grad@x),x0,jac=lambda x:grad,method='SLSQP',bounds=[(-1,1)]*len(K),constraints=[{'type':'ineq','fun':con,'jac':lambda x:con(x,True)}],options={'maxiter':100,'ftol':1e-11,'disp':True})
 K=K0+fit.x*rad;rr=e.rows(allr['t'],allr['u']);vals=e.expression(rr,q);print('KNOT',it,fit.success,fit.fun,'MAX',max(vals),flush=True)
 if max(vals)<-3e-6:break
 ids=np.unique(np.r_[ids,np.argsort(vals)[-3000:]]);x0=fit.x
# Regenerate mesh to include moved kink fibers; keep unscreened if any violation.
base.KNOTS=K;log.ALL_KNOTS=K;log.LOG_KNOTS=K
vv=e.expression(e.grid(dense=True),q);print('MOVED GRID MAX',max(vv),flush=True)
pp=new_to_old(q+1e-8,K);out=dict(dimension=9,parameters=pp.tolist(),v_knots=K.tolist(),log_knots=K.tolist(),monotone_parameters=(q+1e-8).tolist(),certified=False,sampled_max=float(max(vv)))
(D/'knot_fit9.json').write_text(json.dumps(out,indent=2))
