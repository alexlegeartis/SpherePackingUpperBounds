import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import sys,json,math
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog,minimize_scalar
from scipy.special import gammaln,eval_genlaguerre
from flint import fmpq as Q,fmpq_poly as Poly,arb,ctx
from log_hinge import LogHingeExperiment
ctx.prec=256;sys.set_int_max_str_digits(0)
D=Path(__file__).parent;n=int(sys.argv[1]);e=LogHingeExperiment(n);par=np.array(json.loads((D/f'fit{n}.json').read_text())['parameters']);degree=63;s=n/8;alpha=Q(n,2)
rows=np.load(D/f'dim{n}_radial_rows.npz');rows={k:rows[k] for k in rows.files};rows['sign']=np.r_[np.linspace(8,40,2001),rows['sign']]
raw=json.loads((D/f'dim{n}_radial_comparison.json').read_text())['outcomes'][-1]
X=Poly([0,1]);gs=[Poly([1]),alpha-2*X]
for k in range(1,degree):gs.append(((2*k+alpha-2*X)*gs[-1]-(k+alpha-1)*gs[-2])/(k+1))
gs=[g/g(0) for g in gs]
def rat(v):return Q(str(F(float(v))))
def polys(raw):
 p=sum((rat(c)*g for c,g in zip(raw['radial'],gs)),Poly());h=sum(((-1)**k*rat(c)*g for k,(c,g) in enumerate(zip(raw['radial'],gs))),Poly());return p,h
def rr(p,lo):
 zs=p.numer().complex_roots();assert sum(m for z,m in zs)==p.degree();out=[]
 for z,m in zs:
  if z.imag.contains(0) and not z.real<arb(lo):
   assert z.imag.is_exact() and z.imag==0;out.append(float(z.real.mid()))
 return sorted(out)
norm=np.exp(gammaln(np.arange(degree+1)+n/2)-gammaln(n/2)-gammaln(np.arange(degree+1)+1));sg=(-1.)**np.arange(degree+1)
def lag(x):return np.array([eval_genlaguerre(k,n/2-1,2*x)/norm[k] for k in range(degree+1)]).T
def pair(x):
 xx=np.atleast_1d(x);h=np.exp(-e.basis(xx/(np.pi*s))@par);return np.exp(-xx)*(lag(xx)@np.array(raw['radial']))-raw['C']*(2*h-h*h)/e.B
for rnd in range(7):
 p,h=polys(raw)
 for key,pol,lo in [('fourier',h,0),('sign',p,8)]:
  roots=rr(pol,lo);crit=rr(pol.derivative(),lo);pts=crit.copy()
  for a,b in zip(roots[::2],roots[1::2]):pts.extend(np.linspace(max(lo,a-(b-a)),b+(b-a),101))
  rows[key]=np.unique(np.r_[rows[key],pts]);print('ROOTS',rnd,key,len(roots),flush=True)
 xx=np.unique(np.r_[np.linspace(np.pi*s,12,16001),np.pi*s*np.array([1.02,1.05,1.1,1.2,1.4,1.5,1.6,1.65,1.75,2])]);vv=pair(xx);inds=np.where((vv[1:-1]>=vv[:-2])&(vv[1:-1]>=vv[2:]))[0]+1
 for i in inds:
  fit=minimize_scalar(lambda z:-pair(z)[0],bounds=(xx[i-1],xx[i+1]),method='bounded',options={'xatol':1e-14});rows['pair']=np.r_[rows['pair'],fit.x]
 x,y,z=rows['pair'],rows['fourier'],rows['sign'];hh=np.exp(-e.basis(x/(np.pi*s))@par);g=np.exp(x)*(2*hh-hh*hh)/e.B
 A=np.vstack([np.c_[lag(x),-g],np.c_[lag(z),np.zeros(len(z))],np.c_[-lag(y)*sg,np.zeros(len(y))]]);sc=np.maximum(abs(A).max(axis=1),1e-200);A/=sc[:,None]
 fit=linprog(np.ones(degree+2),A_ub=A,b_ub=np.full(len(A),-1e-8),A_eq=np.r_[sg,0][None,:],b_eq=[1],bounds=[(None,None)]*(degree+1)+[(0,None)],method='highs',options={'time_limit':120,'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
 if not fit.success:raise RuntimeError(fit.message)
 raw=dict(radial=fit.x[:-1].tolist(),C=float(fit.x[-1]),objective=float(fit.fun),density_expression=float(fit.fun*np.pi**(n/2)*(s/4)**(n/2)/math.gamma(n/2+1)),certified=False)
 print('LP',rnd,raw['density_expression'],flush=True)
 (D/f'refit{n}.json').write_text(json.dumps(raw,indent=2));np.savez_compressed(D/f'refit_rows{n}.npz',**rows)
 # Make source file consumed by standalone exact export.
 comparison=dict(outcomes=[raw]);(D/f'dim{n}_radial_comparison.json').write_text(json.dumps(comparison,indent=2))
 p,h=polys(raw)
 if not rr(h,0) and not rr(p,8):break
