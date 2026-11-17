# Dominant three-line sign gate, not a standalone density certificate.
import json,math,time
import sympy as s
from mpmath import iv
from pathlib import Path
iv.dps=60;p=Path(__file__).parent
co=json.loads((p/'rank2_candidate_B783489561_20261117.json').read_text())['even_coefficients'];x=s.symbols('x')
Q=s.Poly(sum(s.Rational(str(a))*s.gegenbauer(2*j,23,x)/s.gegenbauer(2*j,23,1)for j,a in enumerate(co))-s.Rational(2,10000),x)
assert Q.eval(0)+s.Rational(1,10000)<0 and s.Poly(Q.as_expr()+s.Rational(1,10000),x).count_roots(0,s.Rational(1,2))==0
assert Q.diff().eval(s.Rational(1,2))>0 and Q.diff().count_roots(s.Rational(1,2),1)==0
cs=[iv.mpf(str(Q.nth(2*j).p))/iv.mpf(str(Q.nth(2*j).q))for j in range(8)]
def up(z):return math.nextafter(float(z.b),math.inf)
def lo(z):return math.nextafter(float(z.a),-math.inf)
def pe(z,der=0):
 mid=iv.mpf((lo(z)+up(z))/2);d=z-mid
 aa=cs if der==0 else [j*cs[j] for j in range(1,8)]
 # Taylor coefficients computed at a point to avoid interval Horner cancellation.
 tt=[sum(aa[j]*math.comb(j,k)*mid**(j-k)for j in range(k,len(aa)))for k in range(len(aa))]
 out=tt[-1]
 for c in tt[-2::-1]:out=out*d+c
 return out
def parts(u,v,gradient=False):
 t=1+u-2*v;a=v*v/u;b=(1-v)**2/t;c=(u-v)**2/(u*t)
 A,B,C=pe(a),pe(b),pe(c);r1=(t/u)**96*iv.exp(-120*(1-2*v));r2=t**96*iv.exp(-120*(t-1));F=A+B*r1+C*r2
 if not gradient:return F
 da,db,dc=pe(a,1),pe(b,1),pe(c,1)
 au=-a/u;av=2*v/u;bu=-b/t;bv=-2*(1-v)/t+2*b/t;cu=c*(2/(u-v)-1/u-1/t);cv=c*(-2/(u-v)+2/t)
 fu=da*au+db*bu*r1+B*r1*96*(1/t-1/u)+dc*cu*r2+C*r2*(96/t-120)
 fv=da*av+db*bv*r1+B*r1*(-192/t+240)+dc*cv*r2+C*r2*(-192/t+240)
 return F,fu,fv
def bound(box):
 ul,ur,vl,vr=box;u=iv.mpf([ul,ur]);v=iv.mpf([vl,vr]);f,fu,fv=parts(u,v,True)
 center=parts((iv.mpf(ul)+iv.mpf(ur))/2,(iv.mpf(vl)+iv.mpf(vr))/2);rad=abs(fu)*(iv.mpf(ur)-ul)/2+abs(fv)*(iv.mpf(vr)-vl)/2
 return min(up(f),up(center+rad))
stack=[(1.,4.,0.,.5)];count=accepted=0;maximum=-math.inf;start=time.time()
while stack and count<100000:
 box=stack.pop();count+=1;b=bound(box)
 if b < -1e-6:accepted+=1;maximum=max(maximum,b)
 else:
  ul,ur,vl,vr=box
  if (ur-ul)/6>vr-vl:m=(ul+ur)/2;stack.extend([(ul,m,vl,vr),(m,ur,vl,vr)])
  else:m=(vl+vr)/2;stack.extend([(ul,ur,vl,m),(ul,ur,m,vr)])
 if count%1000==0:print(count,len(stack),accepted,time.time()-start,box,b,flush=True)
print('DOMINANT GATE',not stack,'nodes',count,'pending',len(stack),'maxupper',maximum,'seconds',time.time()-start,flush=True)

assert not stack, "Incomplete compact sign gate"
