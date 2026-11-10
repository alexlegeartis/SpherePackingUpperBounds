import sys,json,math
sys.set_int_max_str_digits(100000)

from fractions import Fraction as F
from flint import fmpq_poly,fmpq,arb,arb_poly,ctx
from pathlib import Path
ctx.prec=384;path=Path(__file__).resolve().parent;d=json.loads((path/'certificate.json').read_text());src=json.loads((path/d['source']).read_text());assert (src['n'],src['K'],src['L'],src['D'])==(19,5,16,39)
r2=F(19,8);R2=F(57,20);assert F(str(src['ratio']))*r2==R2;ff=lambda x:fmpq(x.numerator,x.denominator)
As=[[[F(z) for z in row] for row in A] for A in d['matrices']];assert len(As)==17
pivots=[]
for l,A in enumerate(As):
 n=6 if l==0 else 5;assert len(A)==n and all(len(row)==n for row in A) and all(A[i][j]==A[j][i] for i in range(n) for j in range(n));S=[row[:] for row in A]
 for k in range(n):
  pivot=S[k][k];assert pivot>0;pivots.append(pivot)
  for i in range(k+1,n):
   for j in range(k+1,n):S[i][j]-=S[i][k]*S[k][j]/pivot
def choose(a,k):return F(math.prod(a-i for i in range(k)),math.factorial(k))
ps=list(map(F,d['physical_polynomial']));hs=list(map(F,d['fourier_polynomial']));qs=[sum(ps[k]*math.factorial(k)*(-1)**j*choose(F(17,2)+k,k-j)/math.factorial(j) for k in range(j,40)) for j in range(40)];assert qs==hs and hs[0]==F(d['normalization']) and hs[0]>0
P=fmpq_poly(list(map(ff,ps)));H=fmpq_poly(list(map(ff,hs)));X=fmpq_poly([0,1])
def chain(p):
 s=[p,p.derivative()]
 while s[-1].degree()>0:
  rem=-(s[-2]%s[-1]);assert rem
  num=rem.numer();s.append(fmpq_poly(num/abs(num.content())))
 return s
def var(s,t=None):
 vals=[g[g.degree()] if t is None else g(ff(t)) for g in s];sg=[int(v>0)-int(v<0) for v in vals if v];return sum(a!=b for a,b in zip(sg,sg[1:]))
# independently generated Machin bounds
def ab(x):
 lo=sum((-1)**j*x**(2*j+1)/(2*j+1) for j in range(40));return lo,lo+x**81/81
al,au=ab(F(1,5));bl,bu=ab(F(1,239));pl=16*al-4*bu;pu=16*au-4*bl
for q in [H,-P(X+ff(R2*pl))]:
 s=chain(q);assert var(s,F(0))==var(s) and q[0]>0
# Explicit normalized Gegenbauer polynomials, not source recurrence.
def rise(a,k):return math.prod(a+i for i in range(k))
G=[]
for l in range(17):
 co=[F(0)]*(l+1);norm=F(rise(17,l),math.factorial(l))
 for j in range(l//2+1):co[l-2*j]=F((-1)**j*rise(F(17,2),l-j)*2**(l-2*j),math.factorial(j)*math.factorial(l-2*j))/norm
 G.append(fmpq_poly(list(map(ff,co))))
nodes=[r2+(R2-r2)*i/32 for i in range(33)]
def w(u):
 q=(u-r2)*4/(R2-r2);i=min(3,int(q));a=[F(0)]*5;a[i]=1-(q-i);a[i+1]=q-i;return a
def cup(u,v):
 z=(u+v-r2)**2/(4*u*v);den=10**40;q=z*den*den;k=math.isqrt(q.numerator//q.denominator)
 while F(k*k,den*den)<z:k+=1
 return F(k,den)
polys={}
for i,u in enumerate(nodes):
 for j in range(i,33):
  a=w(u);b=w(nodes[j]);cs=[]
  for l,A in enumerate(As):
   off=int(l==0);cs.append(sum(a[s]*A[s+off][t+off]*b[t] for s in range(5) for t in range(5)))
  q=sum((ff(c)*g for c,g in zip(cs,G)),fmpq_poly([]));cut=cup(nodes[min(i+1,32)],nodes[min(j+1,32)]);s=chain(q);assert var(s,F(-1))==var(s,cut) and q(-1)<0 and q(ff(cut))<0
  polys[i,j]=True
# Piecewise quadratic coupled kernel, independently reconstruct sparse weights.
AP=arb_poly([arb(ff(x)) for x in ps]);pi=arb.pi();stack=[];cellq=[]
for i in range(4):
 a=r2+(R2-r2)*i/4;b=r2+(R2-r2)*(i+1)/4;ws={i:(ff(b)-X)/ff(b-a),i+1:(X-ff(a))/ff(b-a)}
 q=fmpq_poly([])
 for j,wj in ws.items():
  q+=2*ff(As[0][0][j+1])*wj
  for k,wk in ws.items():q+=ff(As[0][j+1][k+1]+sum(A[j][k] for A in As[1:]))*wj*wk
 cellq.append(arb_poly([arb(c) for c in q.coeffs()]));stack.append((a,b,i,0))
# Direct centered power-series enclosure of the complete physical-plus-kernel function.
passed=0;maxdepth=0
while stack:
 a,b,i,depth=stack.pop();mid=(a+b)/2;rad=(b-a)/2;u=arb(ff(mid)).union(arb(ff(mid)))
 # second-order Taylor and polynomial bounds with independently built Arb polynomials
 xp=pi*arb(ff(mid));h=arb(ff(rad));xx=arb(ff(a)).union(arb(ff(b)))*pi
 P1=AP.derivative()-AP;P2=AP.derivative().derivative()-2*AP.derivative()+AP;q=cellq[i]
 center=AP(xp)*(-xp).exp()+q(u);first=pi*P1(xp)*(-xp).exp()+q.derivative()(u)
 co=P2(arb_poly([xp,pi*h]));mag=sum((abs(c) for c in co.coeffs()),arb(0));second=pi*pi*mag*(-pi*arb(ff(a))).exp()+abs(q.derivative().derivative()[0])
 upper=center+abs(first)*h+second*h*h/2
 if upper<0:passed+=1;maxdepth=max(maxdepth,depth);continue
 assert depth<25;stack.extend([(a,mid,i,depth+1),(mid,b,i,depth+1)])
obj=(ps[0]+As[0][0][0])/hs[0];# Gamma(21/2)=20! sqrt(pi)/(4^10 10!), so density is rational times sqrt(r2).
density_square=obj**2*pu**18*r2**19*4*math.factorial(10)**2/math.factorial(20)**2;assert density_square<F('0.009618351')**2
density=arb(ff(obj))*arb.pi()**9*arb(ff(r2))**F(19,2).__float__()*2*math.factorial(10)/math.factorial(20)
out={'pass':True,'PSD_blocks':17,'positive_pivots':len(pivots),'Fourier_coefficients_match':True,'angular_node_pairs':len(polys),'annular_leaves':passed,'annular_maxdepth':maxdepth,'objective_exact':str(obj),'density_squared_upper_exact':str(density_square),'density_upper_float':float(density.mid()),'strict_density_upper':'0.009618351'}
(path/'independent19_audit.json').write_text(json.dumps(out,indent=2));print({k:v for k,v in out.items() if 'exact' not in k})
