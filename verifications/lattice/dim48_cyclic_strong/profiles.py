if not __debug__:raise RuntimeError('Assertions required')
import sys,json,math,hashlib
from fractions import Fraction as F
from pathlib import Path
import sympy as sp
sys.set_int_max_str_digits(0)
src=Path(sys.argv[1]);raw=src.read_bytes();D=json.loads(raw);n=D['n'];E=list(map(F,D['E']));g=list(map(F,D['g']));J=E.index(F(2));M=len(g);nn=M-J;cuts=D['angular_cuts'];nu=list(map(F,D['nu']));assert len(cuts)==len(nu) and all(v>=0 for v in nu)
x=sp.symbols('x')
def upper(coeff,a,b,parts=32):
 degree=len(coeff)-1;bounds=[]
 for piece in range(parts):
  lo=a+(b-a)*F(piece,parts);width=(b-a)/parts;shift=[sum(coeff[k]*math.comb(k,j)*lo**(k-j) for k in range(j,degree+1)) for j in range(degree+1)]
  bounds.extend(sum(shift[j]*width**j*F(math.comb(i,j),math.comb(degree,j)) for j in range(i+1)) for i in range(degree+1))
 return max(bounds)
for cut in cuts:
 prefix=cut['source_prefix'];U=E[prefix];co=list(map(F,cut['gegenbauer_even']));f=list(map(F,cut['f_y']));assert co[0]==1 and min(co)>=0
 full=sp.Poly(sum(sp.Rational(v.numerator,v.denominator)*sp.gegenbauer(2*k,sp.Rational(n-2,2),x)/sp.gegenbauer(2*k,sp.Rational(n-2,2),1) for k,v in enumerate(co)),x)
 assert all(F(full.nth(2*k))==v for k,v in enumerate(f)) and full.degree()<=2*(len(f)-1)
 assert upper(f,F(0),(U-1)**2)<=0
 f1=sum(f);a0=F(cut['constant']);assert f1==F(cut['f1']) and a0>0 and len(cut['caps'])==J
 for j,(lo,hi) in enumerate(zip(E[:J],E[1:J+1])):
  cmin=max(F(0),(2-hi)/(2*U));cmax=min(1-F(1,2)/U,(2*U-lo)/2);major=max(F(0),upper(f,cmin*cmin,cmax*cmax));row=cut['caps'][j];assert major<=F(row['M'])
  radius=U-lo/4;cos=1-F(1,2)/radius;eta=(U-1)**2/(lo*(4-lo));proj=(cos+eta)/(1-eta);moment=F(1);choices=[]
  for r in range(9):
   th=F(2*r+1,n-1+2*r)
   if proj*proj<th:choices.append((1-proj*proj)/(moment*(th-proj*proj)))
   moment*=th
  cap=2*min(v.numerator//v.denominator for v in choices);assert cap==row['C']
  weight=(2-f1/a0 if j<prefix else 0)-cap*F(row['M'])/a0;assert weight==F(cut['weights'][j])
 print('AA cut PASS',prefix,flush=True)
def B(t):c=1-2/t;return int(n*(1-c*c)/(1-n*c*c))
def A(a,t):
 if a+t<=4:return 0
 c=1-2/t;e=(a+t-4)**2/(4*a*t);cp=(c+e)/(1-e)
 return min(B(t)-1,int((n-1)*(1-cp*cp)/(1-(n-1)*cp*cp))) if (n-1)*cp*cp<1 else B(t)-1
assert len(D['profile_duals'])==J+1
budgets=[]
for idx,entry in enumerate(D['profile_duals']):
 j=idx-1;assert entry['j']==j
 caps=[B(t) if j<0 else A(E[j+1],t) for t in E[J+1:]];caps=[min(caps[k:]) for k in range(nn)]
 rhs=list(map(F,caps));cols=[[F(int(k<=h)) for h in range(nn)] for k in range(nn)]
 for factor in ([F(1),F(21,20)] if j>=0 else [F(1)]):
  for h,t in enumerate(E[J+1:]):
   alpha=factor*(n-1)*(t-2)/(4*(t-1));condition=(j>=0 and (factor>1 or E[j]>=F(191,96)))
   rhs.append(n-1+alpha*(E[j+1]-2) if condition else F(n))
   for k in range(nn):cols[k].append(1-alpha*(E[J+k+1]-2) if k<=h else F(0))
 y=list(map(F,entry['y']));assert len(y)==len(rhs) and min(y)>=0
 for k in range(nn):assert sum(a*b for a,b in zip(y,cols[k]))>=g[J+k]
 bound=sum(a*b for a,b in zip(y,rhs))+(g[j] if j>=0 else F(0))-(sum(v*F(c['weights'][j]) for v,c in zip(nu,cuts)) if j>=0 else 0)
 assert bound==F(entry['bound']);budgets.append(bound)
assert budgets==list(map(F,D['budgets'])) and max(budgets)==F(D['kappa'])
extra=2*sum(v*F(c['constant']) for v,c in zip(nu,cuts));assert extra==F(D['extra_objective'])
ratio=(F(D['p'][0])+F(D['origin_C'])+2*(2**n-1)*max(budgets)+extra)/F(D['q'][0]);density=(F(D['tau'])/8)**(n//2)*ratio/math.factorial(n//2)
assert density==F(D['density']) and density<=F(D['safe_density'])<F('0.0000004084357036')
report=dict(status='PASS',sha256=hashlib.sha256(raw).hexdigest(),AA_cuts=len(cuts),profile_branches=len(budgets),normalization='even exact rational',safe=D['safe_density']);print(report,flush=True);Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2))
