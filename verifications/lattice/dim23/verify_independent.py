"""Independent exact replay of global short-class certificate (generic dimension).
No output claim is accepted as evidence; constraints are reconstructed.
"""
from fractions import Fraction as F
from pathlib import Path
from math import comb,factorial,prod
import json,sympy as S,sys,time,hashlib
if not __debug__:raise RuntimeError('Assertions required')
if hasattr(sys, 'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
assert len(sys.argv)>1, 'Pass the immutable certificate path'
path=Path(sys.argv[1])
raw=path.read_bytes();D=json.loads(raw);n=D['n'];assert n>=4 and F(D['gap'])==0 and D['exp_terms']==128
E=list(map(F,D['E']));J=E.index(F(2));U=E[1];assert E[0]==1 and 1<U<2 and E[-1]==F(9,4) and all(a<b for a,b in zip(E,E[1:]))
x=S.symbols('x');y=S.symbols('y')
def sympoly(values,var=x):return S.Poly(sum(S.Rational(v.numerator,v.denominator)*var**i for i,v in enumerate(values)),var)
def nonnegative_no_roots(poly,a,b):
 # Strict values at an interior point and absence of roots prove strict sign throughout.
 if poly.is_zero:return
 a=S.Rational(a.numerator,a.denominator);b=S.Rational(b.numerator,b.denominator)
 assert poly.eval(a)>=0 and poly.eval(b)>=0
 assert poly.eval((a+b)/2)>0
 # remove possible endpoint factors before checking open-interval roots
 for z in [a,b]:
  while poly.eval(z)==0:poly=poly.exquo(S.Poly(poly.gen-z,poly.gen))
 assert poly.count_roots(a,b)==0
cuts=[]
for rec in D.get('angular_cuts',[]):
 interval=list(map(F,rec['source_interval']));assert rec['n']==n and interval[0]==1
 U=interval[1];assert U in E[1:J] and U<F(5,4)
 if 'source_prefix' in rec:assert rec['source_prefix']==E.index(U)
 co=list(map(F,rec['gegenbauer_even']));assert co[0]==1 and min(co)>=0
 P=S.Poly(sum(S.Rational(c.numerator,c.denominator)*S.gegenbauer(2*k,S.Rational(n-2,2),x)/S.gegenbauer(2*k,S.Rational(n-2,2),1) for k,c in enumerate(co)),x)
 fy=list(map(F,rec['f_y']));assert P==S.Poly(sympoly(fy,y).as_expr().subs(y,x*x),x)
 H=F(P.eval(1));N0=F(rec['constant']);assert H==F(rec['f1']) and N0>0
 nonnegative_no_roots(-sympoly(fy,y),F(0),(U-1)**2)
 assert len(rec['caps'])==J and len(rec['weights'])==J
 weights=[]
 for j,(a,b) in enumerate(zip(E[:J],E[1:J+1])):
  z=rec['caps'][j];assert z['j']==j
  lower=max(F(0),(2-b)/(2*U));upper=min(1-F(1,2)/U,(2*U-a)/2)
  assert list(map(F,z['cosine_squared_interval']))==[lower*lower,upper*upper]
  Mj=F(z['M']);assert Mj>=0
  nonnegative_no_roots(S.Poly(Mj,y)-sympoly(fy,y),lower*lower,upper*upper)
  R2=U-a/4;cs=1-F(1,2)/R2;eps=(U-1)**2/(a*(4-a));c=(cs+eps)/(1-eps);dd=n-1
  assert 0<=cs<1 and eps<1 and R2<1
  choices=[]
  maxr=rec.get('cap_max_r',2);assert isinstance(maxr,int) and 0<=maxr<=8
  for r in range(maxr+1):
   moment=F(prod(2*k+1 for k in range(r)),prod(dd+2*k for k in range(r)))
   threshold=F(2*r+1,dd+2*r)
   if c*c<threshold:choices.append((1-c*c)/(moment*(threshold-c*c)))
  assert choices
  C=2*min(v.numerator//v.denominator for v in choices);assert z['C']>=C
  weight=(2-H/N0 if E[j+1]<=U else F(0))-F(z['C'])*Mj/N0
  assert weight==F(rec['weights'][j]);weights.append(weight)
 cuts.append((N0,weights));print('Independent angular+cap PASS',rec['index'],flush=True)
p=list(map(F,D['p']));q=list(map(F,D['q']));d=list(map(F,D['d']));K=len(p)-1;assert len(q)==len(d)==len(p)
for j in range(K+1):
 rise=prod(F(n,2)+i for i in range(j));den=factorial(j)*rise
 assert p[j]==sum((d[k]*(-1)**j*F(factorial(k),factorial(k-j))/den for k in range(j,K+1)),F(0))
 assert q[j]==sum((d[k]*(-1)**(j+k)*F(factorial(k),factorial(k-j))/den for k in range(j,K+1)),F(0))
tau=F(D['tau']);assert tau>0
Q=sympoly(q);P=sympoly(p)
assert Q.eval(0)>0
if '--bernstein-fourier' in sys.argv:
 from flint import fmpq
 # (1-t)^K q(t/(1-t)) has Bernstein coefficients q_j/binomial(K,j).
 bb=[fmpq(v.numerator,v.denominator*comb(K,j)) for j,v in enumerate(q)]
 pending=[(bb,0)];nodes=0;deep=0
 while pending:
  vv,depth=pending.pop();nodes+=1;deep=max(deep,depth)
  if min(vv)>=0:continue
  assert depth<100 and vv[0]>=0 and vv[-1]>=0
  ll=[vv[0]];rr=[vv[-1]]
  while len(vv)>1:vv=[(a+b)/2 for a,b in zip(vv,vv[1:])];ll.append(vv[0]);rr.append(vv[-1])
  pending.extend([(ll,depth+1),(rr[::-1],depth+1)])
 print('Independent Fourier Mobius Bernstein PASS',nodes,deep,flush=True)
else:
 assert Q.count_roots(0,S.oo)==0
 print('Independent Fourier global sign PASS',flush=True)
assert P.eval(tau*E[-1])<0
if '--bernstein-tail' in sys.argv:
 try:
  from flint import fmpq
 except ImportError:
  fmpq=F
 start=tau*E[-1];end=F(400);T=-P
 while not all(v>=0 for v in T.shift(end).all_coeffs()):end*=2
 Q=T.shift(S.Rational(start.numerator,start.denominator));width=end-start;ww=fmpq(width.numerator,width.denominator)
 aa=[fmpq(int(Q.nth(j).p),int(Q.nth(j).q))*ww**j for j in range(K+1)]
 bb=[sum((aa[j]*fmpq(comb(i,j),comb(K,j)) for j in range(i+1)),fmpq(0)) for i in range(K+1)]
 pending=[(bb,0)];nodes=0;deep=0
 while pending:
  vv,depth=pending.pop();nodes+=1;deep=max(deep,depth)
  if min(vv)>=0:continue
  assert depth<60 and vv[0]>0 and vv[-1]>0
  ll=[vv[0]];rr=[vv[-1]]
  while len(vv)>1:vv=[(a+b)/2 for a,b in zip(vv,vv[1:])];ll.append(vv[0]);rr.append(vv[-1])
  pending.extend([(ll,depth+1),(rr[::-1],depth+1)])
 print('Independent radial tail Bernstein PASS',str(start),str(end),nodes,deep,flush=True)
else:
 assert P.count_roots(tau*E[-1],S.oo)==0
 print('Independent radial tail sign PASS',flush=True)
def shift(poly,a):return [sum((poly[k]*comb(k,j)*a**(k-j) for k in range(j,len(poly))),F(0)) for j in range(len(poly))]
def bern(poly,width):
 degree=len(poly)-1
 return [sum((poly[j]*width**j*F(comb(i,j),comb(degree,j)) for j in range(i+1)),F(0)) for i in range(degree+1)]
def positive_bern(b):
 stack=[(b,0)];depthmax=0
 while stack:
  c,dep=stack.pop();depthmax=max(depthmax,dep)
  if min(c)>=0:continue
  assert dep<24 and c[0]>=0 and c[-1]>=0
  left=[c[0]];right=[c[-1]]
  while len(c)>1:c=[(a+b)/2 for a,b in zip(c,c[1:])];left.append(c[0]);right.append(c[-1])
  stack.extend([(left,dep+1),(right[::-1],dep+1)])
 return depthmax
Gam=list(map(F,D['Gam']));assert len(Gam)==len(E)-1 and min(Gam)>=0
for j in range(len(Gam)):
 a=tau*E[j];width=tau*(E[j+1]-E[j]);poly=[-v for v in shift(p,a)];poly[0]+=Gam[j];poly[1]+=Gam[j]/2;poly[2]+=Gam[j]/8
 positive_bern(bern(poly,width))
print('Independent all radial intervals PASS',len(Gam),flush=True)
def explower(t):return sum((t**k/F(factorial(k)) for k in range(129)),F(0))
g=[v/explower(tau*a/2) for v,a in zip(Gam,E)]
for j in range(len(g)-2,J-1,-1):g[j]=max(g[j],g[j+1])
def B(t):
 c=1-2/t;assert n*c*c<1
 return (n*(1-c*c))//(1-n*c*c)
def A(a,t):
 if a+t<=4:return 0
 c=1-2/t;ca=(a+t-4)**2/(4*a*t);cc=(c+ca)/(1-ca)
 cap=B(t)-1
 if (n-1)*cc*cc<1:cap=min(cap,((n-1)*(1-cc*cc))//(1-(n-1)*cc*cc))
 return cap
# Reconstruct conditioned profile duals independently via suffix sums.
assert n==23 and F(D['short_threshold'])==F(191,96)
from itertools import combinations
excess=[e-2 for e in E[J+1:]];alphas=[F(n-1)*v/(4*(1+v)) for v in excess]
for delta,al in zip(excess,alphas):
 Nmax=F(4*n)*(1+delta)/(4+4*delta-(n-1)*delta*delta)
 for count in range(n+1,Nmax.numerator//Nmax.denominator+1):
  lim=F(count-n)/al
  # Unconditioned rank polynomial is convex; strict negative at 0,
  # nonpositive at lim suffices to exclude every S<lim.
  def unc(Sv):return n*(4*count+4*Sv+F(count,2)*delta*Sv+Sv*Sv/2)-(2*count+Sv)**2
  assert unc(F(0))<0 and unc(lim)<=0 and F(n,2)-1>0
  edges=[(F(0),F(1),F(0)),(F(0),F(1),F(1,96)),(F(1),F(2-count),F(0)),(F(1),F(1),(count-1)*delta),(F(1),F(0),lim)]
  vertices=[]
  for (a,b,c),(d,e,f) in combinations(edges,2):
   det=a*e-b*d
   if not det:continue
   S0=(c*e-b*f)/det;z=(a*f-c*d)/det
   if 0<=z<=F(1,96) and (count-2)*z<=S0<=(count-1)*delta-z and S0<=lim:
    Q=(delta+z)*S0+2*z*z-(count-2)*delta*z
    val=n*(4*count+4*S0+F(count,2)*Q+S0*S0/2)-(2*count+S0)**2
    assert val<=0 and (S0==lim or val<0);vertices.append((S0,z))
  assert vertices
  for a,b,c in edges:
   assert (F(n,2)-1)*b*b-F(n*count,2)*a*b+n*count*a*a>=0
print('Independent unconditional and conditioned weighted rank proofs PASS',flush=True)
num=len(excess);assert len(D['profile_duals'])==J+1
budgets=[]
for j,rec in zip(range(-1,J),D['profile_duals']):
 assert rec['j']==j
 y=list(map(F,rec['y']));assert len(y)==2*num and min(y)>=0
 y0=y[:num];y1=y[num:]
 for k,delta in enumerate(excess):
  coefficient=sum(y0[k:])+sum(v*(1-al*delta) for v,al in zip(y1[k:],alphas[k:]))
  assert coefficient>=g[J+k]
 caps=[B(t) if j<0 else A(E[j+1],t) for t in E[J+1:]]
 for k in range(num-2,-1,-1):caps[k]=min(caps[k],caps[k+1])
 rhs=[F(n-1)+al*(E[j+1]-2) if j>=0 and E[j]>=F(191,96) else F(n) for al in alphas]
 budget=(g[j] if j>=0 else F(0))+sum(v*c for v,c in zip(y0,caps))+sum(v*c for v,c in zip(y1,rhs))
 assert budget==F(rec['bound']);budgets.append(budget)
assert budgets==list(map(F,D['budgets'])) and max(budgets)==F(D['kappa'])
ratio=(p[0]+2*(2**n-1)*max(budgets))/q[0];assert ratio>0
print('Independent profile-dual replay PASS',len(budgets),flush=True)
bound=F(D['density']);assert bound>0
if n%2:
 m=n//2;factor=F(4**(m+1)*factorial(m+1),factorial(2*m+2));sqpi=(tau/8)**n*factor*factor*ratio*ratio
 pilo=16*sum((F((-1)**k,(2*k+1)*5**(2*k+1)) for k in range(4)),F(0))-F(4,239)
 assert pilo>3 and F(D['pi_lower'])==pilo and F(D['density_squared_times_pi'])==sqpi
 assert bound*bound*pilo>=sqpi
else:
 exactdensity=(tau/8)**(n//2)/factorial(n//2)*ratio
 assert bound>=exactdensity
priors={19:F('.00782447361'),17:F('.0138839313'),26:F('.000915119176'),32:F('.000105764483')}
if n in priors:assert bound<priors[n]
assert F(D['safe_density'])>=bound
out={'passed':True,'n':n,'density':str(bound),'improvement_fraction':str(1-bound/priors[n]) if n in priors else None,'sha256':hashlib.sha256(raw).hexdigest(),'angular_cuts':len(cuts),'radial_intervals':len(Gam),'gap':0,'low_density_branch_required':False}
Path(__file__).with_name('audit_conditioned_n'+str(n)+'.json').write_text(json.dumps(out,indent=2));print('ALL INDEPENDENT CHECKS PASS',out,flush=True)
