import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import sys,json,math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import experiment as base
from log_hinge import LogHingeExperiment,V_KNOTS,LOG_KNOTS,sampled_audit
from flint import fmpq as Q,fmpq_poly as P
D=Path(__file__).parent
n=int(sys.argv[1]);j=json.loads((D/f'n{n}_a8_b16.json').read_text())
c=list(map(lambda v:Q(str(F(v))),j['coefficients'])); guard=Q(1,10000);c[0]-=guard;c=[v/c[0] for v in c]
x=P([0,1]);gs=[P([1]),x]
for k in range(1,len(c)-1):gs.append(((2*k+n-2)*x*gs[-1]-k*gs[-2])/(k+n-2))
poly=sum((a*g for a,g in zip(c,gs)),P())
angular={'dimension':n,'a_exact':'8','b_exact':'16','gegenbauer_coefficients':list(map(str,c)),'B_exact':str(sum(c)),'polynomial_coefficients':list(map(str,poly.coeffs()))}
(D/'inputs'/f'angular_dim{n}.json').write_text(json.dumps(angular,indent=2))
e=LogHingeExperiment(n)
rows=e.grid(dense=True)
scale=np.maximum(1e-3,abs(rows['q0'])+abs(rows['q1'])+abs(rows['q2']))
edges=np.r_[1.,V_KNOTS];slopes=8+16*((edges[:-1]-1)+(edges[1:]-1));slopes=np.r_[slopes,40]
p0=np.array(json.loads((D/f'fit{n}.json').read_text())['parameters'])
v=base.VPRICE;bp=e.basis(v);weights=np.exp(-((v-1)/.2)**2);weights/=weights.sum();actrl=2*np.exp(-e.control_phi(v))-np.exp(-2*e.control_phi(v))
def fun(p):
 h=np.exp(-bp@p);return -float(np.sum(weights*(2*h-h*h)/actrl))
def jac(p):
 h=np.exp(-bp@p);return np.sum((weights*2*h*(1-h)/actrl)[:,None]*bp,axis=0)
con={'type':'ineq','fun':lambda p:(-e.expression(rows,p)-1e-5)/scale,'jac':lambda p:-e.expression(rows,p,True)[1]/scale[:,None]}
fit=minimize(fun,p0,jac=jac,method='SLSQP',bounds=[(0,100)]*14,constraints=[con,{'type':'ineq','fun':lambda p:p[:7].sum()-7}],options={'maxiter':120,'ftol':2e-12,'disp':True})
out=dict(dimension=n,parameters=fit.x.tolist(),success=bool(fit.success),certified=False,train_max=float(e.expression(rows,fit.x).max()))
(D/f'fit{n}.json').write_text(json.dumps(out,indent=2));print('PROFILE',out,flush=True)
a=sampled_audit(e,fit.x);(D/f'audit{n}.json').write_text(json.dumps(a,indent=2));print('AUDIT',a,flush=True)
e.radial([dict(index='log',parameters=fit.x.tolist())])
