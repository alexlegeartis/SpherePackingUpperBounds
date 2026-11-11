"""192-bit outward Arb checks of numerical local maxima and exact Gram boundary points."""
if not __debug__:
    raise RuntimeError("Certificate verification requires assertions; do not run Python with -O or -OO.")

from pathlib import Path
from fractions import Fraction
import json,hashlib
from flint import fmpq,arb,ctx
ctx.prec=192
ROOT=Path(__file__).resolve().parent;SRC=ROOT.parent/'density_kissing16_reweight_frozen_20261111'
DATA=json.loads((SRC/'factors.json').read_text(),parse_float=str);W=json.loads((SRC/'candidate.json').read_text(),parse_float=str)['record']
q=lambda x:fmpq(str(Fraction(str(x))))
a=lambda x:arb(q(x))
TH=list(map(a,W['theta']));RR=list(map(a,W['local_angular_coefficients']));FF=[]
for block in DATA['blocks']:
 n=block['size'];L=[[arb(0) for _ in range(n)] for _ in range(n)]
 for i,j,s in block['L_entries']:L[i][j]=a(s)
 F=[[sum((L[k][i]*L[k][j] for k in range(n)),arb(0))+(a(block['lambda_decimal']) if i==j else 0) for j in range(n)] for i in range(n)];FF.append(F)
BB=FF[16];FF=FF[17:]
def gegen(x,n,d):
 p=[arb(1),x]
 for k in range(1,d):p.append(((2*k+n-2)*x*p[-1]-k*p[-2])/(k+n-2))
 return p[:d+1]
def kernel(u,v,t):
 g=t-u*v;z=(1-u*u)*(1-v*v);Q=[arb(1),g]
 for k in range(1,16):Q.append(((2*k+13)*g*Q[-1]-k*z*Q[-2])/(k+13))
 val=BB[1][1]
 for k,F in enumerate(FF):
  if not TH[k]:continue
  pu=gegen(u,16+2*k,16-k);pv=gegen(v,16+2*k,16-k)
  val+=TH[k]*Q[k]*sum((pu[i]*F[i][j]*pv[j] for i in range(len(pu)) for j in range(len(pv))),arb(0))
 return val

def cycle(x):
 u,v,t=x;val=kernel(u,v,t)+kernel(u,t,v)+kernel(v,t,u)
 for z in x:
  p=gegen(z,16,24);val+=sum((RR[i-1]*p[i] for i in range(1,25)),arb(0))
 return val

def det(x):
 u,v,t=x;return 1+2*u*v*t-u*u-v*v-t*t
pts=[('interior_high',['0.18403018856546585','0.18405988347694224','0.18402197292619207']),('rounding_invalid',['-0.9202415586238045','-0.7992488031068492','0.5003']),('other_interior',['0.18994652549234142','-0.23645627418311882','-0.053502856410066246'])]
out=[]
for name,pt in pts:
 x=list(map(a,pt));dd=det(x);val=cycle(x);out.append({'name':name,'exact_decimal_point':pt,'gram_determinant':str(dd),'exact_gram_feasible':bool(dd>=0),'cyclic_enclosure':str(val),'strictly_negative':bool(val<0)})
u=a('-0.9202415586238045');t=arb(fmpq(5003,10000));v=u*t-((1-u*u)*(1-t*t)).sqrt();val=cycle([u,v,t]);out.append({'name':'projected_exact_gram_boundary','u_exact':'-0.9202415586238045','t_exact':'5003/10000','v_definition':'u*t-sqrt((1-u^2)*(1-t^2))','v_enclosure':str(v),'gram_feasible_by_algebraic_identity':True,'coordinates_in_cap':bool(u>=-1 and u<=t and v>=-1 and v<=t),'cyclic_enclosure':str(val),'strictly_negative':bool(val<0)})
res={'arithmetic':'Arb 192 bits','scope':'Point checks only. No continuum certificate. The exactly Gram-invalid optimizer point is not used as proof acceptance.','points':out};(ROOT/'arb_points.json').write_text(json.dumps(res,indent=2));print(json.dumps(res,indent=2))
