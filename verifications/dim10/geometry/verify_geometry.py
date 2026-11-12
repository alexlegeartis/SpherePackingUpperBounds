"""Exact compact cyclic covering for the frozen piecewise-log profile.
At hinge points the AD derivative is an interval containing both one-sided
derivatives. Each hinge is locally Lipschitz; the coordinatewise fundamental
theorem of calculus therefore justifies the same mean-value enclosure.
"""
if not __debug__:raise RuntimeError('Assertions required')
from pathlib import Path
from fractions import Fraction as F
import json,time,hashlib
from flint import arb,arb_poly,fmpq,ctx
ctx.prec=256
D=Path(__file__).resolve().parent
import verify_structure as s
def rat(x):
 q=F(str(x));return fmpq(q.numerator,q.denominator)
AQ=arb_poly([arb(x) for x in s.Q]);DQ=AQ.derivative()
params=list(map(rat,s.fit['parameters']))
VK=list(map(rat,['1.02','1.05','1.1','1.2','1.4','2']))
LK=list(map(rat,['1.2','1.4','1.5','1.6','1.65','1.75']))
def prange(poly,x):
 lo=max(arb(-1),x.lower());hi=min(arb(1),x.upper());assert lo<=hi
 mid=(lo+hi)/2;rad=(hi-lo)/2;co=poly(arb_poly([mid,rad]))
 err=sum((abs(z) for z in co.coeffs()[1:]),arb(0))
 return co[0]+arb(0,err.upper())
class AD:
 def __init__(self,v,dt=0,du=0):self.v=arb(v);self.dt=arb(dt);self.du=arb(du)
 @staticmethod
 def make(x):return x if isinstance(x,AD) else AD(x)
 def __add__(self,o):o=AD.make(o);return AD(self.v+o.v,self.dt+o.dt,self.du+o.du)
 __radd__=__add__
 def __neg__(self):return AD(-self.v,-self.dt,-self.du)
 def __sub__(self,o):return self+-AD.make(o)
 def __rsub__(self,o):return AD.make(o)+-self
 def __mul__(self,o):o=AD.make(o);return AD(self.v*o.v,self.dt*o.v+self.v*o.dt,self.du*o.v+self.v*o.du)
 __rmul__=__mul__
 def inv(self):
  iv=1/self.v;return AD(iv,-self.dt*iv*iv,-self.du*iv*iv)
 def __truediv__(self,o):return self*AD.make(o).inv()
 def sqrt(self):v=self.v.sqrt();return AD(v,self.dt/(2*v),self.du/(2*v))
 def exp(self):v=self.v.exp();return AD(v,v*self.dt,v*self.du)
 def log(self):assert self.v>0;return AD(self.v.log(),self.dt/self.v,self.du/self.v)
 def pos(self):
  if self.v<=0:return AD(0)
  if self.v>=0:return self
  mult=arb(0).union(arb(1))
  return AD(arb(0).union(self.v.upper()),mult*self.dt,mult*self.du)
 def angular(self):d=prange(DQ,self.v);return AD(prange(AQ,self.v),d*self.dt,d*self.du)
def phi(v):
 ans=params[0]*(v-1)+params[7]*v.log()
 for i,k in enumerate(VK):ans+=params[1+i]*(v-k).pos()
 for i,k in enumerate(LK):ans+=params[8+i]*(v/k).log().pos()
 return ans
def evaluate(T,U):
 T=AD(T,1,0);U=AD(U,0,1)
 Y=U*U;Z=1+U*(U-1)+U*(1-2*T)
 for a in [Y,Z]:
  if not a.v>=1:a.v=arb(1).union(a.v.upper())
 W=Z.sqrt();c1=(1-U*T)/W;c2=(U-T)/W
 py=phi(Y);pz=phi(Z)
 return T.angular()+(py-pz).exp()*c1.angular()+(-pz).exp()*c2.angular()
start=time.monotonic();stack=[(fmpq(-1),fmpq(1,2),fmpq(1),fmpq(4),0)]
accept=[];skip=[];tested=0;area=fmpq(0);maxdepth=0
try:
 while stack:
  l,r,v,w,depth=stack.pop();tested+=1
  assert time.monotonic()-start<600 and tested<500000,('BUDGET',tested,len(stack))
  if 2*l*v>1:skip.append([str(x)for x in [l,r,v,w]]);area+=(r-l)*(w-v);continue
  f=evaluate(arb(l).union(arb(r)),arb(v).union(arb(w)))
  assert all(z.is_finite() for z in [f.v,f.dt,f.du]),('NONFINITE',l,r,v,w)
  upper=f.v;method='direct'
  if not upper<0:
   center=evaluate(arb((l+r)/2),arb((v+w)/2)).v
   upper=center+abs(f.dt)*arb((r-l)/2)+abs(f.du)*arb((w-v)/2);method='Lipschitz_mean_value'
  if upper<0:
   accept.append(dict(box=list(map(str,[l,r,v,w])),method=method,upper=str(upper)))
   area+=(r-l)*(w-v);maxdepth=max(maxdepth,depth);continue
  # A positive corner inside the sorted region is a decisive exact failure.
  if depth>=20:
   for tt,uu in [(l,v),(l,w),(r,v),(r,w)]:
    if 2*tt*uu<=1:
     point=evaluate(arb(tt),arb(uu)).v
     assert not point>0,('POSITIVE_CORNER',str(tt),str(uu),str(point))
  assert depth<48,('DEPTH',l,r,v,w)
  if 3*(r-l)>=2*(w-v):
   m=(l+r)/2;stack.extend([(l,m,v,w,depth+1),(m,r,v,w,depth+1)])
  else:
   m=(v+w)/2;stack.extend([(l,r,v,m,depth+1),(l,r,m,w,depth+1)])
 assert area==fmpq(9,2)
 out=dict(status='PASS',scope='Angular kernel geometry only, not radial or density proof',accepted=len(accept),discarded=len(skip),tested=tested,maxdepth=maxdepth,exact_area=str(area),seconds=time.monotonic()-start,
  input_hashes={p:hashlib.sha256((D/p).read_bytes()).hexdigest() for p in ['fit.json','angular.json']})
 (D/'partition.json').write_text(json.dumps(dict(accepted=accept,discarded=skip,area=str(area)))+'\n')
except Exception as e:
 out=dict(status='INCOMPLETE',reason=str(e),tested=tested,accepted=len(accept),remaining=len(stack),seconds=time.monotonic()-start)
 (D/'unresolved.json').write_text(json.dumps(dict(current=list(map(str,[l,r,v,w,depth])),stack=[list(map(str,z)) for z in stack]))+'\n')
 (D/'RESULT_GEOMETRY.json').write_text(json.dumps(out,indent=2)+'\n');print(out,flush=True);raise
(D/'RESULT_GEOMETRY.json').write_text(json.dumps(out,indent=2)+'\n');print(out,flush=True)
