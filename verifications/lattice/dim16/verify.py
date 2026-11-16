if not __debug__:
 raise RuntimeError("Verification requires assertions enabled; do not use python -O")
import json,math,time,sys
from fractions import Fraction as F
import sympy as s
from pathlib import Path
path=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('certificate.json')
c=json.load(open(path));assert c['n']==16 and F(c['tau0'])==16 and c['exp_terms']==128
p=list(map(F,c['p']));q=list(map(F,c['q']));d=list(map(F,c['d']));E=list(map(F,c['E']));tau=F(16);n=16;J=E.index(F(2));M=len(E)-1;N=len(p)-1;mult=2*(2**n-1)
assert len(q)==len(p)==len(d) and E[0]==1 and E[-1]==F(9,4) and all(a<b for a,b in zip(E,E[1:]))
for j in range(N+1):
 assert p[j]==sum((d[k]*(-1)**j*F(math.comb(k+7,k-j),math.factorial(j)*math.comb(k+7,k)) for k in range(j,N+1)),F(0))
 assert q[j]==sum((d[k]*(-1)**(j+k)*F(math.comb(k+7,k-j),math.factorial(j)*math.comb(k+7,k)) for k in range(j,N+1)),F(0))
x=s.Symbol('x');P=s.Poly.from_list([s.Rational(v.numerator,v.denominator) for v in p[::-1]],x);Q=s.Poly.from_list([s.Rational(v.numerator,v.denominator) for v in q[::-1]],x)
assert Q.eval(0)>0 and Q.count_roots(0,s.oo)==0
assert P.eval(36)<0 and P.count_roots(36,s.oo)==0
print('Exact Fourier/tail signs PASS',flush=True)
assert F(16)*F(80,81)/F(65,81)<20
for count in [17,18,19]:
 S=F(4*(count-16),3);C=4*count-F(count*count,4)+(4-F(count,8))*S+F(7,16)*S*S
 assert C<0 and 4-F(count,8)>0
print('Weighted-class lemma constants PASS',flush=True)
def exp_lower(x):
 term=s=F(1)
 for i in range(1,129):term*=x/i;s+=term
 return s
def shift(p,a):
 c=[p[-1]]
 for z in p[-2::-1]:
  new=[a*c[0]+z]+[c[i-1]+a*c[i] for i in range(1,len(c))]+[c[-1]];c=new
 return c
def bern(c,w):
 c=[v*w**i for i,v in enumerate(c)]
 return [sum((c[j]*F(math.comb(i,j),math.comb(N,j)) for j in range(i+1)),F(0)) for i in range(N+1)]
def check(b,depth=0):
 if min(b)>=0:return True,depth
 if b[0]<0 or b[-1]<0:return False,depth
 if depth>=16:return False,depth
 left=[b[0]];right=[b[-1]]
 for r in range(N):b=[(a+c)/2 for a,c in zip(b,b[1:])];left.append(b[0]);right.append(b[-1])
 u,d=check(left,depth+1)
 if not u:return u,d
 v,e=check(right[::-1],depth+1);return v,max(d,e)
gap=F(c['gap']);alpha=F(c['alpha']);assert gap==F(4,15) and alpha==F(4,5)
for count in [17,18,19]:
 S=gap+F(count-16)/alpha
 assert 4*count-F(count*count,4)+(4-F(count,8))*S+F(7,16)*S*S<0
Gam=list(map(F,c['Gam']));Mweight=F(c['Mweight']);assert len(Gam)==M and min(Gam)>=0 and Mweight>0
for b in range(M):
 a=tau*E[b];w=tau*(E[b+1]-E[b]);sp=shift(p,a);poly=[-v for v in sp];poly[0]+=Gam[b];poly[1]+=Gam[b]/2;poly[2]+=Gam[b]/8
 assert check(bern(poly,w))[0],('G interval',b)
 if b>=J:
  fac=Mweight*exp_lower(a/2);l0=1+2*alpha-alpha*a/tau;l1=-alpha/tau
  poly=[-v for v in sp];poly[0]+=fac*l0;poly[1]+=fac*(l0/2+l1);poly[2]+=fac*(l0/8+l1/2);poly[3]+=fac*l1/8
  assert check(bern(poly,w))[0],('weighted interval',b)
print('All spatial intervals PASS',flush=True)
g=[ga/exp_lower(tau*e/2) for ga,e in zip(Gam,E)]
for b in range(M-2,J-1,-1):g[b]=max(g[b],g[b+1])
g.append(F(0))
def B(s):c=1-2/s;return (n*(1-c*c))//(1-n*c*c)
def A(a,s):
 if a+s<=4:return 0
 c=1-2/s;ca2=(a+s-4)**2/(4*a*s);cp=(c+ca2)/(1-ca2)
 return min(B(s)-1,((n-1)*(1-cp*cp))//(1-(n-1)*cp*cp)) if (n-1)*cp*cp<1 else B(s)-1
nn=M-J
ae=[E[b+1]-2 for b in range(J,M)]
caps=[[1+A(E[j+1],E[b+1]) for b in range(J,M)] for j in range(J)]
for row in caps:
 for b in range(nn-2,-1,-1):row[b]=min(row[b],row[b+1])
duals={}
for z in c['profile_duals']:
 key=(z['j'],z['h']);assert key not in duals;duals[key]=z['y']
expected=set();bounds=[];affected=[]
def diffs(seq):return [seq[0]]+[b-a for a,b in zip(seq,seq[1:])]
for j in range(J):
 sa=E[j+1]-2;first=diffs([min(v-1,15) for v in caps[j]])
 ftrace=sa+sum((v*a for v,a in zip(first,ae)),F(0))
 if caps[j][-1]<16 or ftrace>=gap:
  bounds.append(g[j]+sum(((g[b]-g[b+1])*(caps[j][b-J]-1) for b in range(J,M)),F(0)));continue
 affected.append(j)
 light=diffs([min(v-1,14) for v in caps[j]])
 bound=g[j]+sum((v*w for v,w in zip(light,g[J:M])),F(0))
 for h in range(nn):
  rhs3=15*ae[h]+sa-gap
  if caps[j][h]<16 or rhs3<0:continue
  expected.add((j,h));assert (j,h) in duals
  y=[F(0)]*(nn+3);seen=set()
  for pos,value in duals[(j,h)]:
   assert isinstance(pos,int) and 0<=pos<nn+3 and pos not in seen;seen.add(pos);y[pos]=F(value);assert y[pos]>=0
  total=F(0);lhs=[F(0)]*nn
  for b in range(nn-1,-1,-1):total+=y[b];lhs[b]=total
  for b in range(nn):
   if b<h:lhs[b]+=y[nn]+y[nn+2]*(ae[h]-ae[b])
   if b<=h:lhs[b]-=y[nn+1]
   assert lhs[b]>=g[J+b],('dual feasibility',j,h,b)
  ub=g[j]+sum(((caps[j][b]-1)*y[b] for b in range(nn)),F(0))+15*y[nn]-15*y[nn+1]+rhs3*y[nn+2]
  bound=max(bound,ub)
 bounds.append(bound)
assert set(duals)==expected and affected==c['affected']
print('All',len(expected),'exact profile duals PASS',flush=True)
empty=(16-alpha*gap)*Mweight;k=max([empty]+bounds)
density=(tau/8)**8/F(math.factorial(8))*(p[0]+mult*k)/q[0]
assert empty==F(c['empty_bound']) and k==F(c['kappa']) and density==F(c['density'])
assert 0<density<=F('0.018655117428')
# Low-density alternative rho<=4096 is below the displayed bound using pi<22/7.
assert F(22,7)**8/F(16*math.factorial(8))<F('0.018655117428')
print('PASS: density <=0.018655117428, using the proved near-frame dichotomy')
