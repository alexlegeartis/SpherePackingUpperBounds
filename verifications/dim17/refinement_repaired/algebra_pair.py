import sys,json,math
sys.set_int_max_str_digits(100000)

from fractions import Fraction as F
from flint import fmpq_poly,fmpq,arb,arb_poly,ctx
from pathlib import Path
assert __debug__
ctx.prec=384;path=Path(__file__).resolve().parent;d=json.loads((path/'certificate.json').read_text());assert (d['dimension'],d['radial_nodes'],d['angular_degree'],d['radial_degree'])==(17,7,16,31)
r2=F(d['exclusion_squared']);R2=F(d['truncation_squared']);assert r2==F(17,8) and R2==F("2.9749999999999996") and r2<R2<2*r2
ff=lambda x:fmpq(x.numerator,x.denominator)
As=[[[F(z) for z in row] for row in A] for A in d['matrices']];assert len(As)==17
pivots=[]
for l,A in enumerate(As):
 n=8 if l==0 else 7;assert len(A)==n and all(len(row)==n for row in A) and all(A[i][j]==A[j][i] for i in range(n) for j in range(n));S=[row[:] for row in A]
 for k in range(n):
  pivot=S[k][k];assert pivot>0;pivots.append(pivot)
  for i in range(k+1,n):
   for j in range(k+1,n):S[i][j]-=S[i][k]*S[k][j]/pivot
def choose(a,k):return F(math.prod(a-i for i in range(k)),math.factorial(k))
ps=list(map(F,d['physical_polynomial']));hs=list(map(F,d['fourier_polynomial']));qs=[sum(ps[k]*math.factorial(k)*(-1)**j*choose(F(15,2)+k,k-j)/math.factorial(j) for k in range(j,32)) for j in range(32)];assert qs==hs and hs[0]==F(d['normalization']) and hs[0]>0
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
 lo=sum((-1)**j*x**(2*j+1)/(2*j+1) for j in range(32));return lo,lo+x**65/65
al,au=ab(F(1,5));bl,bu=ab(F(1,239));pl=16*al-4*bu;pu=16*au-4*bl
for q in [H]:
 s=chain(q);assert var(s,F(0))==var(s) and q[0]>0
# Piecewise quadratic coupled kernel, independently reconstruct sparse weights.
AP=arb_poly([arb(ff(x)) for x in ps]);pi=arb.pi();stack=[];cellq=[]
for i in range(6):
 a=r2+(R2-r2)*i/6;b=r2+(R2-r2)*(i+1)/6;ws={i:(ff(b)-X)/ff(b-a),i+1:(X-ff(a))/ff(b-a)}
 q=fmpq_poly([])
 for j,wj in ws.items():
  q+=2*ff(As[0][0][j+1])*wj
  for k,wk in ws.items():q+=sum((ff(A[j+int(ell==0)][k+int(ell==0)])*(X/ff(r2))**ell for ell,A in enumerate(As)),fmpq_poly([]))*wj*wk
 cellq.append(arb_poly([arb(c) for c in q.coeffs()]));stack.append((a,b,i,0))
# Direct centered power-series enclosure of the complete physical-plus-kernel function.
passed=0;maxdepth=0
while stack:
 a,b,i,depth=stack.pop();mid=(a+b)/2;rad=(b-a)/2;u=arb(ff(mid)).union(arb(ff(mid)))
 # second-order Taylor and polynomial bounds with independently built Arb polynomials
 xp=pi*arb(ff(mid));h=arb(ff(rad));xx=arb(ff(a)).union(arb(ff(b)))*pi
 P1=AP.derivative()-AP;P2=AP.derivative().derivative()-2*AP.derivative()+AP;q=cellq[i]
 center=AP(xp)*(-xp).exp()+q(u);first=pi*P1(xp)*(-xp).exp()+q.derivative()(u)
 co=P2(arb_poly([xp,pi*h]));mag=sum((abs(c) for c in co.coeffs()),arb(0));second=pi*pi*mag*(-pi*arb(ff(a))).exp()+abs(q.derivative().derivative()(arb(ff(a)).union(arb(ff(b)))))
 upper=center+abs(first)*h+second*h*h/2
 if upper<0:passed+=1;maxdepth=max(maxdepth,depth);continue
 assert depth<25;stack.extend([(a,mid,i,depth+1),(mid,b,i,depth+1)])

from repair_checks import check_repair
out=check_repair(d,ps,hs,As,r2,R2,pl,pu,passed,maxdepth,len(pivots))
(path/"algebra_pair_audit.json").write_text(json.dumps(out,indent=2));print(out,flush=True)
