import sys,os,json,math,time,hashlib
from pathlib import Path
from fractions import Fraction as F
from flint import arb,arb_poly,ctx
if not __debug__:raise RuntimeError("Assertions must remain enabled")
ctx.prec=1024
def aa(x):
 x=F(x);return arb(x.numerator)/x.denominator
class Jet:
 def __init__(self,v,g=None,h=None):
  self.v=v if isinstance(v,arb)else arb(v);self.g=g if g is not None else[arb(0)]*3;self.h=h if h is not None else[[arb(0)]*3 for _ in range(3)]
 def __add__(self,b):
  b=b if isinstance(b,Jet)else Jet(b);return Jet(self.v+b.v,[x+y for x,y in zip(self.g,b.g)],[[self.h[i][j]+b.h[i][j]for j in range(3)]for i in range(3)])
 __radd__=__add__
 def __neg__(self):return Jet(-self.v,[-x for x in self.g],[[-x for x in row]for row in self.h])
 def __sub__(self,b):return self+-b if isinstance(b,Jet)else self+(-b)
 def __rsub__(self,b):return -self+b
 def __mul__(self,b):
  b=b if isinstance(b,Jet)else Jet(b);return Jet(self.v*b.v,[self.g[i]*b.v+self.v*b.g[i]for i in range(3)],[[self.h[i][j]*b.v+self.v*b.h[i][j]+self.g[i]*b.g[j]+self.g[j]*b.g[i]for j in range(3)]for i in range(3)])
 __rmul__=__mul__
 def __truediv__(self,b):
  b=b if isinstance(b,Jet)else Jet(b);v=1/b.v;iv=Jet(v,[-x*v*v for x in b.g],[[2*b.g[i]*b.g[j]*v**3-b.h[i][j]*v*v for j in range(3)]for i in range(3)]);return self*iv
 def log(self):return Jet(self.v.log(),[x/self.v for x in self.g],[[self.h[i][j]/self.v-self.g[i]*self.g[j]/self.v**2 for j in range(3)]for i in range(3)])
 def exp(self):
  v=self.v.exp();return Jet(v,[v*x for x in self.g],[[v*(self.h[i][j]+self.g[i]*self.g[j])for j in range(3)]for i in range(3)])

out=Path(sys.argv[1]).resolve();r=json.loads((out/'angular.json').read_text());cs=[aa(F(x)) for x in r['polynomial_coefficients']]
assert F(r['profile_a_exact'])==-192 and F(r['profile_b_exact'])==216
assert ctx.prec>=384
assert len(cs)==r['degree']+1 and cs[-1]!=0
ap=arb_poly(cs)
def pe2(unused,z):
 shifted=ap(arb_poly([z.v.mid(),arb(1)]));delta=z.v-z.v.mid()
 value=shifted(delta);der=shifted.derivative()(delta);sec=shifted.derivative().derivative()(delta)
 return Jet(value,[der*x for x in z.g],[[sec*z.g[i]*z.g[j]+der*z.h[i][j]for j in range(3)]for i in range(3)])
def parts(x,u,mode):
 t=x/(2*u) if mode else x
 gap=1-x if mode else 1-2*u*x
 w2=u*u+gap;w=(w2.log()/2).exp()
 r1=(192*(1+gap/(u*u)).log()-216*gap).exp()
 r2=(192*w2.log()-216*(w2-1)).exp()
 return pe2(cs,t)+r1*pe2(cs,(1-u*t)/w)+r2*pe2(cs,(u-t)/w)
start=time.time();count=passed=0;worst=arb('-inf');maxdep=0
stack=[(mode,(F(0) if mode else F(-1),F(1) if mode else F(0),F(1),F(8)),0) for mode in [0,1]]
while stack:
 mode,b,d=stack.pop();count+=1;m=[(b[2*i]+b[2*i+1])/2 for i in range(2)];rad=[(b[2*i+1]-b[2*i])/2 for i in range(2)]
 assert all(aa(z).is_exact() for z in m+rad)
 args=[Jet(arb(aa(m[i]),aa(rad[i])),[arb(int(i==j))for j in range(3)])for i in range(2)]
 f=parts(*args,mode);c=parts(*[Jet(aa(m[i]),[arb(int(i==j))for j in range(3)])for i in range(2)],mode)
 bound=c.v+sum((abs(c.g[i])*aa(rad[i])for i in range(2)),arb(0))+sum((abs(f.h[i][j])*aa(rad[i]*rad[j]/2)for i in range(2)for j in range(2)),arb(0))
 first=c.v+sum((abs(f.g[i])*aa(rad[i])for i in range(2)),arb(0))
 if f.v<0 or bound<0 or first<0:passed+=1;maxdep=max(maxdep,d);continue
 if c.v>0:
  print('EXACT POSITIVE',mode,list(map(str,m)),str(c.v),flush=True);raise SystemExit(2)
 assert d<80,(mode,b)
 scores=[float((abs(c.g[i])*aa(rad[i])+sum((abs(f.h[i][j])*aa(rad[i]*rad[j]) for j in range(2)),arb(0))).upper()) for i in range(2)]
 ax=max(range(2),key=lambda i:scores[i]) if all(math.isfinite(v) for v in scores) else max(range(2),key=lambda i:rad[i]);j=2*ax;mid=(b[j]+b[j+1])/2;l=list(b);rr=list(b);l[j+1]=mid;rr[j]=mid;stack.extend([(mode,tuple(l),d+1),(mode,tuple(rr),d+1)])
 if count%1000==0:print('TAYLOR',count,passed,len(stack),mode,list(map(str,b)),time.time()-start,flush=True)
assert count==2*passed-2
res=dict(angular_sha256=hashlib.sha256((out/'angular.json').read_bytes()).hexdigest(),status='COMPACT TAYLOR PASS',precision=ctx.prec,processed=count,leaves=passed,maxdepth=maxdep,seconds=time.time()-start,domain='u in[1,8]; t in[-1,0] or t=v/(2u), v in[0,1]')
(out/'taylor_fast_angular_result.json').write_text(json.dumps(res,indent=2));print(res,flush=True)
