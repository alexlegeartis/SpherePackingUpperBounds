import json,sys,hashlib,math
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations
import sympy as S
from flint import fmpq as Q,arb,ctx
from interval import upper
sys.set_int_max_str_digits(0);ctx.prec=512
if not __debug__:raise RuntimeError('assertions required')
path=Path(sys.argv[1]);raw=path.read_bytes();D=json.loads(raw);n=D['n'];assert n==46
tau=F(D['tau']);E=list(map(F,D['E']));J=E.index(F(2));M=len(E)-1;assert E[0]==1 and E[-1]<=F(9,4) and all(a<b for a,b in zip(E,E[1:]));p=list(map(F,D['p']));q=list(map(F,D['q']));d=list(map(F,D['d']));K=len(p)-1;assert len(q)==len(d)==len(p)
for j in range(K+1):
 den=math.factorial(j)*math.prod(F(n,2)+i for i in range(j));z=[d[k]*(-1)**j*F(math.factorial(k),math.factorial(k-j))/den for k in range(j,K+1)]
 assert p[j]==sum(z) and q[j]==sum(v*(-1)**k for k,v in zip(range(j,K+1),z))
def shift(p,a):
 out=[p[-1]]
 for v in p[-2::-1]:out=[a*out[0]+v]+[out[i-1]+a*out[i] for i in range(1,len(out))]+[out[-1]]
 return out
def ray(c,a):
 c=shift(c,a);assert c[0]>0;bb=[Q(v.numerator,v.denominator*math.comb(K,j)) for j,v in enumerate(c)];stack=[(bb,0)];nodes=0
 while stack:
  v,h=stack.pop();nodes+=1
  if min(v)>=0:continue
  assert h<100 and v[0]>=0 and v[-1]>=0
  left=[v[0]];right=[v[-1]]
  while len(v)>1:v=[(a+b)/2 for a,b in zip(v,v[1:])];left.append(v[0]);right.append(v[-1])
  stack.extend([(left,h+1),(right[::-1],h+1)])
 return nodes
print('independent global rays',ray(q,F(0)),ray([-v for v in p],tau*E[-1]),flush=True)
cy=D['cyclic_angular'];angular_path=Path('cyclic/angular.json');assert hashlib.sha256(angular_path.read_bytes()).hexdigest()=='c64c10c88bb6af72608674e5bd34ea4b4bf298eeaa7cfe780a0e3dcf189e9bdf';assert cy==json.loads(angular_path.read_bytes());assert cy==json.load(open('cyclic/angular.json'));assert json.load(open('cyclic/meanvalue_angular_result.json'))['pass']
B=F(cy['B_exact']);qm=sum((-1)**k*F(v) for k,v in enumerate(cy['gegenbauer_coefficients']));C=F(D['C']);pb=F(cy['profile_b_exact']);power=-int(F(cy['profile_a_exact']));assert C>=0 and B>0 and 0<1+qm/B<1 and pb>power>0
def ar(f):return arb(f.numerator)/f.denominator
origin_ratio=F(D['origin_ratio']);assert origin_ratio==F(7,10);origin_C=F(D['origin_C']);assert origin_C==C*origin_ratio**2
h4=arb(4)**power*(-ar(pb)*3).exp();assert ar(1+qm/B)+2*h4<2*ar(origin_ratio)
# Source angular cuts: independent coefficient reconstruction and exact root exclusion.
x=S.symbols('x');y=S.symbols('y')
def polynomial(c,var):return S.Poly(sum(S.Rational(v.numerator,v.denominator)*var**j for j,v in enumerate(c)),var)
def nonneg(poly,a,b):
 if poly.is_zero:return
 a=S.Rational(a);b=S.Rational(b);assert poly.eval(a)>=0 and poly.eval(b)>=0 and poly.eval((a+b)/2)>0
 for t in [a,b]:
  while poly.eval(t)==0:poly=poly.exquo(S.Poly(poly.gen-t,poly.gen))
 assert poly.count_roots(a,b)==0
for cut in D['angular_cuts']:
 pref=cut['source_prefix'];U=E[pref];assert cut['n']==n and list(map(F,cut['source_interval']))==[1,U] and U<F(5,4)
 co=list(map(F,cut['gegenbauer_even']));assert co[0]==1 and min(co)>=0
 P=S.Poly(sum(S.Rational(v)*S.gegenbauer(2*k,S.Rational(n-2,2),x)/S.gegenbauer(2*k,S.Rational(n-2,2),1) for k,v in enumerate(co)),x)
 fy=list(map(F,cut['f_y']));assert P==S.Poly(polynomial(fy,y).as_expr().subs(y,x*x),x);N0=F(cut['constant']);H=sum(co);assert N0>0 and H==F(cut['f1']);nonneg(-polynomial(fy,y),F(0),(U-1)**2)
 assert len(cut['caps'])==len(cut['weights'])==J
 for j,z in enumerate(cut['caps']):
  a,b=E[j:j+2];lo=max(F(0),(2-b)/(2*U));hi=min(1-F(1,2)/U,(2*U-a)/2);assert z['j']==j and list(map(F,z['cosine_squared_interval']))==[lo*lo,hi*hi]
  mj=F(z['M']);assert mj>=0;nonneg(S.Poly(mj,y)-polynomial(fy,y),lo*lo,hi*hi)
  R=U-a/4;ep=(U-1)**2/(a*(4-a));c=(1-F(1,2)/R+ep)/(1-ep);dd=n-1;assert R<1 and ep<1;opts=[]
  for r in range(cut.get('cap_max_r',2)+1):
   mu=F(math.prod(2*k+1 for k in range(r)),math.prod(dd+2*k for k in range(r)));t=F(2*r+1,dd+2*r)
   if c*c<t:
    v=(1-c*c)/(mu*(t-c*c));opts.append(2*(v.numerator//v.denominator))
  assert opts and z['C']>=min(opts)
  assert F(cut['weights'][j])==(2-H/N0 if j<pref else 0)-F(z['C'])*mj/N0
print('independent AA cuts PASS',len(D['angular_cuts']),flush=True)
# Complete gap-free saved partitions; own centered-Taylor enclosures.
g=list(map(F,D['g']));assert len(g)==M and min(g)>=0 and all(g[k]>=g[k+1] for k in range(J,M-1));assert len(D['radial_leaves'])==M
nodes=0
for j,keys in enumerate(D['radial_leaves']):
 intervals=sorted((F(int(k,2) if k else 0,2**len(k)),F((int(k,2) if k else 0)+1,2**len(k))) for k in keys);end=F(0)
 for a,b in intervals:assert a==end;end=b
 assert end==1
 stack=[(E[j]+a*(E[j+1]-E[j]),E[j]+b*(E[j+1]-E[j]),0) for a,b in intervals]
 while stack:
  a,b,depth=stack.pop();nodes+=1
  vals=upper(p,tau,C,B,qm,power,pb,a,b,origin_ratio)
  if all(v<=ar(g[j]) for v in vals):continue
  assert depth<18,('interval',j,a,b)
  mid=(a+b)/2;stack.extend([(a,mid,depth+1),(mid,b,depth+1)])
 if j%30==0:print('independent mixed intervals',j,'nodes',nodes,flush=True)
print('independent mixed envelopes PASS',M,nodes,flush=True)
def Bcap(t):
 c=1-2/t;assert n*c*c<1;return int(n*(1-c*c)/(1-n*c*c))
def Acap(a,t):
 if a+t<=4:return 0
 c=1-2/t;e=(a+t-4)**2/(4*a*t);cp=(c+e)/(1-e);v=Bcap(t)-1
 if (n-1)*cp*cp<1:v=min(v,int((n-1)*(1-cp*cp)/(1-(n-1)*cp*cp)))
 return v
exc=[t-2 for t in E[J+1:]];nn=len(exc);alphas=[F(n-1)*t/(4*(1+t)) for t in exc]
# Independently prove all needed weighted prefixes, including all-short extension.
for delta,alpha in zip(exc,alphas):
 nc=4*n*(1+delta)/(4+4*delta-(n-1)*delta*delta)
 for count in range(n+1,int(nc)+1):
  lim=F(count-n)/alpha
  fn=lambda s:n*(4*count+4*s+F(count,2)*delta*s+s*s/2)-(2*count+s)**2
  assert fn(0)<0 and fn(lim)<=0
  for factor,Z,clip in [(F(1),F(1,96),False),(F(21,20),delta,True)]:
   lim=F(count-n)/(factor*alpha);zm=min(Z,delta)
   def possible(z):
    beta=(delta-z)**2/(4*(2-z)*(2+delta));c=(delta/(2+delta)+beta)/(1-beta)
    return c*c>=F(count-n,(n-1)*(count-2))
   if clip:
    if not possible(F(0)):continue
    if not possible(zm):
     lo=F(0);hi=zm
     for _ in range(18):
      mid=(lo+hi)/2
      if possible(mid):lo=mid
      else:hi=mid
     zm=hi
   lines=[(F(0),F(1),F(0)),(F(0),F(1),zm),(F(1),F(2-count),F(0)),(F(1),F(1),(count-1)*delta),(F(1),F(0),lim)];nv=0
   for (a,b,c),(e,f,h) in combinations(lines,2):
    det=a*f-b*e
    if not det:continue
    s=(c*f-b*h)/det;z=(a*h-c*e)/det
    if 0<=z<=zm and (count-2)*z<=s<=(count-1)*delta-z and s<=lim:
     nv+=1;qq=(delta+z)*s+2*z*z-(count-2)*delta*z;v=n*(4*count+4*s+F(count,2)*qq+s*s/2)-(2*count+s)**2;assert v<=0 and (s==lim or v<0)
   assert nv
   for a,b,c in lines:assert (F(n,2)-1)*b*b-F(n*count,2)*a*b+n*count*a*a>=0
print('independent weighted theorems PASS',flush=True)
nu=list(map(F,D['nu']));cuts=D['angular_cuts'];assert len(nu)==len(cuts) and all(v>=0 for v in nu);budgets=[];assert len(D['profile_duals'])==J+1
for j,rec in zip(range(-1,J),D['profile_duals']):
 assert rec['j']==j;yy=list(map(F,rec['y']));groups=2 if j<0 else 3;assert len(yy)==groups*nn and min(yy)>=0
 y0,y1=yy[:nn],yy[nn:2*nn];y2=yy[2*nn:]
 for k,t in enumerate(exc):
  val=sum(y0[k:])+sum(v*(1-al*t) for v,al in zip(y1[k:],alphas[k:]))+sum(v*(1-F(21,20)*al*t) for v,al in zip(y2[k:],alphas[k:]));assert val>=g[J+k]
 caps=[Bcap(t) if j<0 else Acap(E[j+1],t) for t in E[J+1:]]
 for k in range(nn-2,-1,-1):caps[k]=min(caps[k],caps[k+1])
 rhs1=[n-1+al*(E[j+1]-2) if j>=0 and E[j]>=F(191,96) else F(n) for al in alphas];rhs2=[n-1+F(21,20)*al*(E[j+1]-2) for al in alphas] if j>=0 else []
 val=(g[j] if j>=0 else 0)+sum(v*c for v,c in zip(y0,caps))+sum(v*c for v,c in zip(y1,rhs1))+sum(v*c for v,c in zip(y2,rhs2))-(sum(v*F(c['weights'][j]) for v,c in zip(nu,cuts)) if j>=0 else 0)
 assert val==F(rec['bound']);budgets.append(val)
assert budgets==list(map(F,D['budgets'])) and max(budgets)==F(D['kappa']);extra=2*sum(v*F(c['constant']) for v,c in zip(nu,cuts));assert extra==F(D['extra_objective'])
ratio=(p[0]+origin_C+2*(2**n-1)*max(budgets)+extra)/q[0];assert ratio>0
density=(tau/8)**(n//2)/math.factorial(n//2)*ratio;assert density>0 and density==F(D['density']);safe=F(D['safe_density']);assert density<safe
out={'passed':True,'sha256':hashlib.sha256(raw).hexdigest(),'density':str(density),'safe_density':str(safe),'strict':True,'angular_cuts':len(cuts),'radial_bins':M,'mixed_enclosure_nodes':nodes,'class_branches':len(budgets),'source_cyclic_gates':'own46 certified cyclic kernel'}
Path('second_result.json').write_text(json.dumps(out,indent=2));print('FULL SECOND INDEPENDENT HYBRID PASS',out,flush=True)
