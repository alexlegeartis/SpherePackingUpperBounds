import sys,json,math,time

from flint import fmpq,fmpq_poly,arb,ctx
from fractions import Fraction as F
ctx.prec=256
r=json.load(open(sys.argv[1]));n=r['dimension'];K=r['radial_nodes'];r2=F(r['exclusion_squared']);R2=F(r['truncation_squared']);As=[[[F(z)for z in row]for row in A]for A in r['matrices']];P=fmpq_poly(list(map(fmpq,r['physical_polynomial'])));P1=P.derivative()-P;P2=P.derivative().derivative()-2*P.derivative()+P;U=fmpq_poly([0,1]);pi=arb.pi()
def ev(p,x):
 y=arb(0)
 for c in reversed(p.coeffs()):y=y*x+arb(c)
 return y
def mag(p,x,h):
 co=[arb(0)]*(p.degree()+1)
 for c in reversed(p.coeffs()):
  for j in range(p.degree(),0,-1):co[j]=co[j]*x+co[j-1]
  co[0]=co[0]*x+arb(c)
 z=arb(0)
 for c in reversed(co):z=z*h+abs(c)
 return z
stack=[];qs=[]
for i in range(K-1):
 a=r2+(R2-r2)*i/(K-1);b=r2+(R2-r2)*(i+1)/(K-1);ww=[fmpq_poly([])for j in range(K)];ww[i]=(fmpq(str(b))-U)/fmpq(str(b-a));ww[i+1]=(U-fmpq(str(a)))/fmpq(str(b-a));q=fmpq_poly([])
 for j in range(K):q+=2*fmpq(str(As[0][0][j+1]))*ww[j]
 for j in range(K):
  for k in range(K):q+=fmpq(str(As[0][j+1][k+1]+sum(A[j][k]for A in As[1:])))*ww[j]*ww[k]
 qs.append(q);stack.append((a,b,i,0))
checked=0;maxdepth=0;fails=[];t=time.time()
while stack:
 a,b,i,d=stack.pop();mid=(a+b)/2;hh=(b-a)/2;u=arb(fmpq(str(mid)));h=arb(fmpq(str(hh)));xx=pi*u;ex=(-xx).exp();q=qs[i];g=ev(P,xx)*ex+ev(q,u);gp=pi*ev(P1,xx)*ex+ev(q.derivative(),u);M2=pi*pi*mag(P2,xx,pi*h)*(-pi*arb(fmpq(str(a)))).exp()+abs(arb(q.derivative().derivative()[0]));bound=g+abs(gp)*h+M2*h*h/2;checked+=1;maxdepth=max(maxdepth,d)
 if bound<0:continue
 if d>=28 or checked>100000:
  fails.append([str(a),str(b),str(bound),str(g)]);break
 stack.extend([(a,mid,i,d+1),(mid,b,i,d+1)])
print('ANNULAR',len(fails)==0,'boxes',checked,'depth',maxdepth,'seconds',time.time()-t,'failure',fails[:1],flush=True)

assert not fails and not stack
print('Full annular inequality: PASS',flush=True)

