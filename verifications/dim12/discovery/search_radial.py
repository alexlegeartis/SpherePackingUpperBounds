"""Pure-Laguerre curved-profile cyclic radial LP; numerical discovery only.
No source polynomial or fixed early-negativity cutoff is imposed.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,json,math,argparse
from pathlib import Path
from fractions import Fraction as F
from scipy.optimize import linprog
from scipy.special import eval_genlaguerre,gammaln
D=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--dimension',type=int,default=15);ap.add_argument('--degree',type=int,default=47);ap.add_argument('--s',default='15/8');ap.add_argument('--a',default='12');ap.add_argument('--b',default='16');ap.add_argument('--cut',type=float,default=40);ap.add_argument('--B',type=float);ap.add_argument('--angular');ap.add_argument('--margin',type=float,default=1e-9);ap.add_argument('--tag',default='');ap.add_argument('--method',default='highs');a=ap.parse_args();n=a.dimension;d=a.degree;s=float(F(a.s));aa=float(F(a.a));bb=float(F(a.b))
ang=Path(a.angular) if a.angular else D/'angular_exact.json'
j=json.loads(ang.read_text());B=a.B if a.B else float(F(j['B_exact']));name=a.tag or 'radial'
assert a.cut>np.pi*s
norm=np.exp(gammaln(np.arange(d+1)+n/2)-gammaln(n/2)-gammaln(np.arange(d+1)+1));sign=(-1.)**np.arange(d+1)
def lag(X):return np.array([eval_genlaguerre(k,n/2-1,2*X)/norm[k] for k in range(d+1)]).T
# Include a fixed Fourier/ray grid and a much denser near-contact pair grid.
y=np.r_[np.linspace(0,4*d+n+20,6001),np.geomspace(4*d+n+20,2500,801)]
x=np.unique(np.r_[np.linspace(np.pi*s,a.cut,1501),np.pi*s+np.geomspace(1e-8,min(2.,a.cut-np.pi*s),501)])
z=np.r_[np.linspace(a.cut,4*d+n+20,4001),np.geomspace(4*d+n+20,2500,801)]
extra=D/(name+'_cuts.json')
if extra.exists():
 e=json.loads(extra.read_text());y=np.r_[y,e.get('fourier',[])];z=np.r_[z,[v for v in e.get('physical',[]) if v>=a.cut]];x=np.r_[x,[v for v in e.get('pair',[])+e.get('physical',[]) if np.pi*s<=v<a.cut]]
zeta=x/(np.pi*s)-1;h=np.exp(-aa*zeta-bb*zeta*zeta);g=np.exp(x)*(2*h-h*h)/B
A=np.vstack([np.c_[lag(x),-g],np.c_[lag(z),np.zeros(len(z))],np.c_[-lag(y)*sign,np.zeros(len(y))]])
scale=np.maximum(abs(A).max(axis=1),1e-200);A/=scale[:,None]
eq=np.r_[sign,0.];c=np.ones(d+2)
r=linprog(c,A_ub=A,b_ub=np.full(len(A),-a.margin),A_eq=eq[None,:],b_eq=[1.],bounds=[(None,None)]*(d+1)+[(0,None)],method=a.method,options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9,'time_limit':120})
if not r.success:print(r.message,flush=True);raise SystemExit(1)
v=r.x;density=r.fun*np.pi**(n/2)*(s/4)**(n/2)/math.gamma(n/2+1)
out={'dimension':n,'degree':d,'s_exact':a.s,'a_exact':a.a,'b_exact':a.b,'cut_X':a.cut,'B':B,'angular_source':str(ang),'margin':a.margin,'certified':False,'objective':float(r.fun),'density_candidate':float(density),'radial':v[:-1].tolist(),'C':float(v[-1]),'max_scaled_residual':float(max(A@v)),'normalization_residual':float(eq@v-1),'variable_count':len(v),'row_count':len(A)}
(D/(name+'.json')).write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items() if k!='radial'},flush=True)
