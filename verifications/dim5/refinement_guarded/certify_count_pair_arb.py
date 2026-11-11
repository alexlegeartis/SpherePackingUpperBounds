
if not __debug__: raise RuntimeError("Verification must not run under python -O.")
from pathlib import Path
import sys,json,time,math
from fractions import Fraction as F

from flint import fmpq,fmpq_poly,arb,ctx
ctx.prec=256
path=sys.argv[1];data=json.load(open(path));src=json.load(open((Path(path).resolve().parent/data['source']) if not Path(data['source']).is_absolute() else data['source']));R=F(str(src['R']));p=src['p'];K=src['K'];Q=[[[fmpq(z)for z in row]for row in A]for A in data['matrices']];X=fmpq_poly([0,1]);z=(2*X-fmpq(str(R+1)))/fmpq(str(R-1));T=[fmpq_poly([1]),z]
for i in range(1,p):T.append(2*z*T[-1]-T[-2])
q=fmpq_poly([])
for i in range(p+1):q+=2*Q[0][0][i+1]*T[i]
for k in range(K+1):
 o=(k==0)
 for i in range(p+1):
  for j in range(p+1):q+=Q[k][i+o][j+o]*T[i]*T[j]
m=fmpq(src["dimension"],2)
ls=[fmpq_poly([1]),m-2*X]
for k in range(1,31):ls.append(((2*k+m-2*X)*ls[-1]-(k+m-1)*ls[-2])/fmpq(k+1))
P=sum((fmpq(str(F(str(c))))*ls[k]/ls[k][0]for k,c in enumerate(src['x'][:32])),fmpq_poly([]));P1=P.derivative()-P;P2=P.derivative().derivative()-2*P.derivative()+P;q1=q.derivative();q2=q1.derivative();pi=arb.pi()
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
radii=[F(str(z)) for z in src.get('count_radii',[])];weights=[F(str(z)) for z in src.get('count_weights',[])];edges=sorted(set([F(1),R]+[r for r in radii if 1<r<R]));stack=[(a,b,0) for a,b in zip(edges[:-1],edges[1:])];checked=0;maxdepth=0;fail=None;t0=time.time()
while stack:
 a,b,d=stack.pop();m=(a+b)/2;h=(b-a)/2;ma=arb(fmpq(str(m)));ha=arb(fmpq(str(h)));xx=pi*ma*ma;ex=(-xx).exp()
 allowance=sum((w for r,w in zip(radii,weights) if m<r),F(0));g=ev(P,xx)*ex+ev(q,ma)-arb(fmpq(str(allowance)));gp=2*pi*ma*ev(P1,xx)*ex+ev(q1,ma)
 xm=pi*arb(fmpq(str((a*a+b*b)/2)));xh=pi*arb(fmpq(str((b*b-a*a)/2)))
 M2=(2*pi*mag(P1,xm,xh)+4*pi*pi*arb(fmpq(str(b*b)))*mag(P2,xm,xh))*(-pi*arb(fmpq(str(a*a)))).exp()+mag(q2,ma,ha)
 bd=g+abs(gp)*ha+M2*ha*ha/2;checked+=1;maxdepth=max(maxdepth,d)
 if bd<0:continue
 if d>=35 or checked>200000:
  fail={'a':str(a),'b':str(b),'bound':str(bd),'center':str(g)};break
 stack.extend([(a,m,d+1),(m,b,d+1)])
out={'verified':fail is None and not stack,'checked':checked,'max_depth':maxdepth,'seconds':time.time()-t0,'failure':fail,'source':path,'method':'256-bit Arb second-order Taylor bounds on exact rational radial subdivisions'}
print(json.dumps(out),flush=True);json.dump(out,open(path.replace('_psd_exact.json','_count_pair_arb.json'),'w'),indent=2)
