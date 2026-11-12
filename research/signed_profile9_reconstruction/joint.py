import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from scipy.optimize import minimize
from flint import fmpq as Q,fmpq_poly as Poly
import experiment as base
import log_hinge as log
from basis import basis,new_to_old,old_to_new
D=Path(__file__).parent
j=json.loads((D/'fit9.json').read_text());vk=np.array(j['v_knots']);lk=np.array(j['log_knots']);p=np.array(j['parameters'])
K=np.unique(np.r_[vk,lk,[1.003,1.006,1.01,1.015,1.03,1.075,1.15,1.3,1.45,1.55,1.7,1.8,1.9]])
pp=np.r_[p[0],[sum(p[1:1+len(vk)][vk==k]) for k in K],p[1+len(vk)],[sum(p[2+len(vk):][lk==k]) for k in K]]
_,q0=old_to_new(pp,K,K);q0=np.maximum(q0,0)
base.KNOTS=K;log.ALL_KNOTS=K;log.LOG_KNOTS=K
class E(log.LogHingeExperiment):
 @staticmethod
 def basis(v):return basis(v,K)
e=E(9);c0=e.c.copy();L=len(c0);cs=np.maximum(c0[1:],1);qs=np.maximum(q0,1)
def gegen(t):
 t=np.asarray(t);g=[np.ones_like(t),t]
 for k in range(1,L-1):g.append(((2*k+7)*t*g[-1]-k*g[-2])/(k+7))
 return np.stack(g,axis=-1)
raw=json.loads((D/'refit9.json').read_text());z=np.load(D/'duals9.npz');x=z['x'];v=x/(np.pi*9/8);h=np.exp(-basis(v,K)@q0);mult=-z['dual']*raw['C']*np.exp(x)
gq=(mult*2*h*(1-h)/e.B)@basis(v,K);gB=np.sum(mult*(2*h-h*h)/e.B**2);grad=np.r_[gB*cs,gq*qs];grad/=np.linalg.norm(grad)
x0=np.r_[c0[1:]/cs,q0/qs];allr=e.grid(dense=True);vals=e.expression(allr,q0);rng=np.random.default_rng(99);ids=np.unique(np.r_[np.argsort(vals)[-8000:],rng.choice(len(vals),2000,replace=False)]);cap=gegen(np.linspace(-1,.5,4001))
# Tail slope >=7 for v>=16.
tail=np.zeros(len(q0));tail[-1]=1
for i in range(len(q0)):
 unit=np.eye(len(q0))[i];pold=new_to_old(unit,K);tail[i]+=sum(pold[1+len(K):])/16
bound=[(max(0,float(z)-.02),float(z)+.02) for z in x0[:L-1]]+[(max(0,float(z)-2),float(z)+2) for z in x0[L-1:]]
for it in range(7):
 r={k:v[ids] for k,v in allr.items()};t,u=r['t'],r['u'];w=np.sqrt(r['z']);G0=gegen(t);G1=gegen((1-u*t)/w);G2=gegen((u-t)/w);scale=np.maximum(1e-3,abs(r['q0'])+abs(r['q1'])+abs(r['q2']))
 def con(x,jac=False):
  c=np.r_[1,x[:L-1]*cs];q=x[L-1:]*qs;e1=np.exp(-r['diff']@q);e2=np.exp(-r['end']@q);q1=G1@c;q2=G2@c;val=G0@c+e1*q1+e2*q2
  if jac:
   Jc=(G0+e1[:,None]*G1+e2[:,None]*G2)[:,1:]*cs;Jq=(-e1[:,None]*q1[:,None]*r['diff']-e2[:,None]*q2[:,None]*r['end'])*qs
   return np.vstack([-np.c_[Jc,Jq]/scale[:,None],np.c_[-cap[:,1:]*cs,np.zeros((len(cap),len(q))) ]])
  return np.r_[(-val-1e-5)/scale,-cap@c-1e-5]
 fit=minimize(lambda x:float(grad@x),x0,jac=lambda x:grad,method='SLSQP',bounds=bound,constraints=[{'type':'ineq','fun':con,'jac':lambda x:con(x,True)},{'type':'ineq','fun':lambda x:np.r_[tail@(x[L-1:]*qs)-7,x[-1]*qs[-1]-7]}],options={'maxiter':200,'ftol':2e-12,'disp':True})
 c=np.r_[1,fit.x[:L-1]*cs];q=fit.x[L-1:]*qs;e.c=c;e.B=float(sum(c));vals=e.expression(e.rows(allr['t'],allr['u']),q);print('JOINT',it,fit.success,fit.fun,'B',e.B,'MAX',float(vals.max()),flush=True)
 (D/'checkpoint.json').write_text(json.dumps(dict(iteration=it,c=c.tolist(),q=q.tolist(),knots=K.tolist(),maximum=float(vals.max())),indent=2))
 if vals.max()<-3e-6:break
 ids=np.unique(np.r_[ids,np.argsort(vals)[-3000:]]);x0=fit.x
# Positive q repair guards exact monotonic cone, then translate to signed coefficients.
q+=1e-8;pp=new_to_old(q,K);cc=[Q(str(F(float(z)))) for z in c];cc[0]-=Q(1,100000);cc=[z/cc[0] for z in cc]
X=Poly([0,1]);gs=[Poly([1]),X]
for k in range(1,L-1):gs.append(((2*k+7)*X*gs[-1]-k*gs[-2])/(k+7))
qq=sum((z*g for z,g in zip(cc,gs)),Poly());angular=dict(dimension=9,a_exact='8',b_exact='16',gegenbauer_coefficients=list(map(str,cc)),polynomial_coefficients=list(map(str,qq.coeffs())),B_exact=str(sum(cc)))
out=dict(dimension=9,parameters=pp.tolist(),v_knots=K.tolist(),log_knots=K.tolist(),monotone_parameters=q.tolist(),certified=False)
(D/'joint_fit9.json').write_text(json.dumps(out,indent=2));(D/'joint_angular9.json').write_text(json.dumps(angular,indent=2))
