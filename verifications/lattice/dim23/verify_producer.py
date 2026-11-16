"""Standalone exact generic first-n profile certificate replay. Python+SymPy."""
import sys,json,math
from fractions import Fraction as F
if hasattr(sys,"set_int_max_str_digits"):sys.set_int_max_str_digits(0)
import sympy as sp
if not __debug__:raise RuntimeError('Assertions required; do not use python -O')
D=json.load(open(sys.argv[1]));n=D['n'];gap=F(D['gap']);tau=F(D['tau']);E=list(map(F,D['E']));J=E.index(F(2));M=len(E)-1;nn=M-J
p=list(map(F,D['p']));q=list(map(F,D['q']));d=list(map(F,D['d']));N=len(p)-1;alpha=sp.Rational(n,2)-1
assert 3<=n<=48 and E[-1]<=F(9,4)
assert E[0]==1 and E[J]==2 and all(a<b for a,b in zip(E,E[1:]))
for j in range(N+1):
 z=[d[k]*(-1)**j*F(sp.binomial(k+alpha,k-j)/(math.factorial(j)*sp.binomial(k+alpha,k))) for k in range(j,N+1)]
 assert p[j]==sum(z,F(0))
 assert q[j]==sum((v*(-1)**k for k,v in zip(range(j,N+1),z)),F(0))
for name,poly,left in [('q',q,F(0)),('minus_p',[-v for v in p],tau*E[-1])]:
 co=[sum(poly[k]*math.comb(k,j)*left**(k-j) for k in range(j,N+1))/math.comb(N,j) for j in range(N+1)]
 den=math.lcm(*(v.denominator for v in co));stack=[([v.numerator*(den//v.denominator) for v in co],0)]
 while stack:
  rr,depth=stack.pop()
  if min(rr)>0:continue
  assert rr[0]>0 and rr[-1]>0 and depth<100
  ll=[rr[0]<<N];right=[rr[-1]<<N]
  for level in range(1,N+1):
   rr=[aa+bb for aa,bb in zip(rr,rr[1:])];ll.append(rr[0]<<(N-level));right.append(rr[-1]<<(N-level))
  for arr in (ll,right[::-1]):
   gcd=math.gcd(*arr);stack.append(([v//gcd for v in arr],depth+1))
 print(name,'global Mobius Bernstein PASS',flush=True)
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

Gam=list(map(F,D['Gam']));assert len(Gam)==M
for b in range(M):
 a=tau*E[b];width=tau*(E[b+1]-E[b]);co=[-v for v in shift(p,a)];co[0]+=Gam[b];co[1]+=Gam[b]/2;co[2]+=Gam[b]/8
 assert Gam[b]>=0 and check(bern(co,width))[0],('interval',b)
g=[ga/exp_lower(tau*e/2) for ga,e in zip(Gam,E)]
for b in range(M-2,J-1,-1):g[b]=max(g[b],g[b+1])
def B(s):
 c=1-2/s;assert 1-n*c*c>0
 return int(n*(1-c*c)//(1-n*c*c))
def A(a,s):
 if a+s<=4:return 0
 c=1-2/s;ca=(a+s-4)**2/(4*a*s);cp=(c+ca)/(1-ca)
 return min(B(s)-1,int((n-1)*(1-cp*cp)//(1-(n-1)*cp*cp))) if (n-1)*cp*cp<1 else B(s)-1

# Replay the exact conditioned weighted theorem for every used cutoff.
from itertools import combinations
excess=[e-2 for e in E[J+1:]]
for dd in excess:
 al=F(n-1)*dd/(4*(1+dd));nc=F(4*n)*(1+dd)/(4+4*dd-(n-1)*dd*dd)
 for count in range(n+1,int(nc)+1):
  limit=F(count-n)/al
  lines=[(F(0),F(1),F(0)),(F(0),F(1),F(1,96)),(F(1),F(2-count),F(0)),(F(1),F(1),(count-1)*dd),(F(1),F(0),limit)]
  vertices=set()
  for (aa,bb,cc),(ee,ff,gg) in combinations(lines,2):
   det=aa*ff-bb*ee
   if det:
    S=(cc*ff-bb*gg)/det;z=(aa*gg-cc*ee)/det
    if 0<=z<=F(1,96) and (count-2)*z<=S<=(count-1)*dd-z and S<=limit:vertices.add((S,z))
  assert vertices
  for S,z in vertices:
   Qmax=(dd+z)*S+2*z*z-(count-2)*dd*z
   val=n*(4*count+4*S+F(count,2)*Qmax+S*S/2)-(2*count+S)**2
   assert val<=0 and (S==limit or val<0)
  for aa,bb,cc in lines:
   hS,hz=bb,-aa
   assert (F(n,2)-1)*hS*hS+F(n*count,2)*hS*hz+n*count*hz*hz>=0
assert len(D['profile_duals'])==J+1
budgets=[]
for j,dual in zip(range(-1,J),D['profile_duals']):
 assert dual['j']==j
 caps=[B(E[b+1]) if j<0 else A(E[j+1],E[b+1]) for b in range(J,M)]
 for h in range(nn-2,-1,-1):caps[h]=min(caps[h],caps[h+1])
 mat=[[F(int(k<=h)) for k in range(nn)] for h in range(nn)];rhs=list(map(F,caps))
 for h,dd in enumerate(excess):
  al=F(n-1)*dd/(4*(1+dd));mat.append([1-al*e if k<=h else F(0) for k,e in enumerate(excess)])
  rhs.append(F(n-1)+al*(E[j+1]-2) if j>=0 and E[j]>=F(191,96) else F(n))
 y=list(map(F,dual['y']));assert len(y)==2*nn and min(y)>=0
 assert all(sum(y[i]*mat[i][k] for i in range(2*nn))>=g[J+k] for k in range(nn))
 bound=(g[j] if j>=0 else F(0))+sum(v*b for v,b in zip(y,rhs));assert bound==F(dual['bound']);budgets.append(bound)
assert budgets==list(map(F,D['budgets'])) and max(budgets)==F(D['kappa'])
ratio=(p[0]+2*(2**n-1)*max(budgets))/q[0];assert ratio>0 and n%2==1
mh=n//2;factor=F(4**(mh+1)*math.factorial(mh+1),math.factorial(2*mh+2));sqpi=(tau/8)**n*factor**2*ratio**2
pi_lower=16*sum((F((-1)**j,(2*j+1)*5**(2*j+1)) for j in range(4)),F(0))-F(4,239)
assert sqpi==F(D['density_squared_times_pi']) and pi_lower==F(D['pi_lower'])
assert F(D['safe_density'])**2*pi_lower>=sqpi and F(D['safe_density'])>=F(D['density'])
print('FULL CONDITIONED CERTIFICATE PASS',n,D['safe_density'],flush=True)
