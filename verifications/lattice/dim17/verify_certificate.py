"""Standalone exact generic first-n profile certificate replay. Python+SymPy."""
import sys,json,math
from fractions import Fraction as F
import sympy as sp
if not __debug__:raise RuntimeError('Assertions required; do not use python -O')
D=json.load(open(sys.argv[1]));n=D['n'];gap=F(D['gap']);tau=F(D['tau']);E=list(map(F,D['E']));J=E.index(F(2));M=len(E)-1;nn=M-J
p=list(map(F,D['p']));q=list(map(F,D['q']));d=list(map(F,D['d']));N=len(p)-1;alpha=sp.Rational(n,2)-1
assert E[0]==1 and E[J]==2 and all(a<b for a,b in zip(E,E[1:]))
for j in range(N+1):
 z=[d[k]*(-1)**j*F(sp.binomial(k+alpha,k-j)/(math.factorial(j)*sp.binomial(k+alpha,k))) for k in range(j,N+1)]
 assert p[j]==sum(z,F(0))
 assert q[j]==sum((v*(-1)**k for k,v in zip(range(j,N+1),z)),F(0))
x=sp.Symbol('x')
for name,c,left in [('q',q,F(0)),('minus_p',[-v for v in p],tau*E[-1])]:
 poly=sp.Poly.from_list([sp.Rational(v.numerator,v.denominator) for v in c[::-1]],x)
 assert poly.eval(left)>0 and poly.count_roots(left,sp.oo)==0
 print(name,'global sign PASS',flush=True)
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

# All angular checks use exact rational Bernstein upper bounds in y=x^2.
def upper_small(poly,a,b,parts):
 deg=len(poly)-1;out=None
 for h in range(parts):
  l=a+(b-a)*F(h,parts);w=(b-a)/parts;co=shift(poly,l)
  v=max(sum(co[j]*w**j*F(math.comb(i,j),math.comb(deg,j)) for j in range(i+1)) for i in range(deg+1))
  out=v if out is None else max(out,v)
 return out
cuts=D['angular_cuts'];nu=list(map(F,D['nu']));assert len(cuts)==len(nu) and min(nu)>=0
for cut in cuts:
 assert cut['n']==n and list(map(F,cut['source_interval']))==[F(1),E[1]]
 U=E[1];poly=list(map(F,cut['f_y']));co=list(map(F,cut['gegenbauer_even']));assert co[0]==1 and min(co)>=0
 reconstructed=sp.Poly(sum(sp.Rational(v.numerator,v.denominator)*sp.gegenbauer(2*k,sp.Rational(n-2,2),x)/sp.gegenbauer(2*k,sp.Rational(n-2,2),1) for k,v in enumerate(co)),x)
 assert all(F(reconstructed.nth(2*j))==poly[j] for j in range(len(poly))) and reconstructed.degree()==2*(len(poly)-1)
 parts=cut['bernstein_parts'];assert upper_small(poly,F(0),(U-1)**2,parts)<=0
 f1=sum(poly);assert f1==F(cut['f1']);N0=F(cut['constant']);assert N0>0
 assert len(cut['caps'])==J and len(cut['weights'])==J
 for j,(a,b,cap) in enumerate(zip(E[:J],E[1:J+1],cut['caps'])):
  lo=max(F(0),(2-b)/(2*U));hi=min(1-F(1,2)/U,(2*U-a)/2)
  assert cap['j']==j and list(map(F,cap['cosine_squared_interval']))==[lo*lo,hi*hi]
  Mj=F(cap['M']);assert Mj>=0 and Mj>=upper_small(poly,lo*lo,hi*hi,parts)
  R2=U-a/4;c=1-F(1,2)/R2;eta=(U-1)**2/(a*(4-a));cp=(c+eta)/(1-eta);dim=n-1;options=[]
  if dim*cp*cp<1:options.append(dim*(1-cp*cp)/(1-dim*cp*cp))
  if (dim+2)*cp*cp<3:options.append(dim*(dim+2)*(1-cp*cp)/(3-(dim+2)*cp*cp))
  assert options;Cj=2*min(v.numerator//v.denominator for v in options);assert Cj==cap['C']
  assert F(cut['weights'][j])==(2-f1/N0 if j==0 else 0)-Cj*Mj/N0
 print('angular cut PASS',cut['index'],flush=True)
assert gap==0
budgets=[];basebudgets=[]
for j in range(-1,J):
 caps=[B(E[b+1]) if j<0 else A(E[j+1],E[b+1]) for b in range(J,M)]
 for h in range(nn-2,-1,-1):caps[h]=min(caps[h],caps[h+1])
 base=F(0) if j<0 else g[j];last=0
 for cap,w in zip(caps,g[J:]):base+=(cap-last)*w;last=cap
 reward=F(0) if j<0 else sum(v*F(c['weights'][j]) for v,c in zip(nu,cuts))
 basebudgets.append(base);budgets.append(base-reward)
assert basebudgets==list(map(F,D['basebudgets'])) and budgets==list(map(F,D['budgets']))
kappa=max(budgets);assert kappa==F(D['kappa']);extra=2*sum(v*F(c['constant']) for v,c in zip(nu,cuts));assert extra==F(D['extra_objective'])
ratio=(p[0]+2*(2**n-1)*kappa+extra)/q[0];assert ratio>0
if n%2==0:
 density=(tau/8)**(n//2)/math.factorial(n//2)*ratio;assert density==F(D['density'])
else:
 mhalf=n//2;factor=F(4**(mhalf+1)*math.factorial(mhalf+1),math.factorial(2*mhalf+2));squared_pi=(tau/8)**n*factor**2*ratio**2
 pi_lower=16*sum((F((-1)**j,(2*j+1)*5**(2*j+1)) for j in range(4)),F(0))-F(4,239)
 assert pi_lower>3 and F(D['pi_lower'])==pi_lower and squared_pi==F(D['density_squared_times_pi'])
 density=F(D['density']);assert density>0 and density*density*pi_lower>=squared_pi
print('PASS SCHUR dimension',n,'intervals',M,'angular cuts',len(cuts),'density',float(density))
