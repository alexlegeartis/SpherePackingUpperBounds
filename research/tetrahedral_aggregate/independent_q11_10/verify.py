from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import json,math,time
from flint import fmpq as Q,fmpq_poly as Poly,arb,ctx
if not __debug__:raise RuntimeError('Assertions required')
ctx.prec=160;D=Path(__file__).resolve().parent;s=json.loads((D/'certificate.json').read_text());cs=list(map(Q,s['coefficients']));assert s['dimension']==22 and Q(s['q'])==Q(11,10) and cs[0]==1 and min(cs)>=0
# Independent explicit Gegenbauer binomial expansion, lambda=(21-2)/2.
x=Poly([0,1]);lam=Q(19,2)
def rising(a,k):return math.prod((a+j for j in range(k)),start=Q(1))
G=[]
for k in range(len(cs)):
 C=sum(((-1)**j*rising(lam,k-j)*(2*x)**(k-2*j)/(math.factorial(j)*math.factorial(k-2*j))for j in range(k//2+1)),Poly([]));G.append(C/C(1))
P=sum((c*g for c,g in zip(cs,G)),Poly([]));assert P==Poly(list(map(Q,s['polynomial']))) and P(1)==sum(cs)==Q(s['B_exact']);assert P(1)<Q(508698083741,100000000)
# Exact Sturm far sign; no stored success flag used.
a=Q(-1);b=Q(6,17);assert b==1-Q(11,10)/(2*(Q(11,10)-Q(1,4)));assert P(a)<0 and P(b)<0
seq=[P,P.derivative()]
while seq[-1].degree()>0:
 r=-(seq[-2]%seq[-1]);assert r;seq.append(r)
def V(t):
 signs=[int(p(t)>0)-int(p(t)<0)for p in seq if p(t)];return sum(a!=b for a,b in zip(signs,signs[1:]))
assert V(a)==V(b)
class A:
 def __init__(self,v,d=None):self.v=v if isinstance(v,arb)else arb(v);self.d=d or [arb(0)]*6
 def __add__(self,o):
  o=o if isinstance(o,A)else A(o);return A(self.v+o.v,[a+b for a,b in zip(self.d,o.d)])
 __radd__=__add__
 def __neg__(self):return A(-self.v,[-d for d in self.d])
 def __sub__(self,o):return self+-o
 def __rsub__(self,o):return -self+o
 def __mul__(self,o):
  o=o if isinstance(o,A)else A(o);return A(self.v*o.v,[a*o.v+self.v*b for a,b in zip(self.d,o.d)])
 __rmul__=__mul__
 def inv(self):return A(1/self.v,[-d/self.v**2 for d in self.d])
 def __truediv__(self,o):return self*(o.inv()if isinstance(o,A)else A(o).inv())
 def sqrt(self):
  v=self.v.sqrt();return A(v,[d/(2*v)for d in self.d])
edges=list(combinations(range(4),2));idx={e:i for i,e in enumerate(edges)}
def value(box,gradient):
 vals=[]
 for i,(lo,hi)in enumerate(box):
  ds=[arb(int(i==j))for j in range(6)]if gradient else [arb(0)]*6
  vals.append(A(arb((lo+hi)/2,(hi-lo)/2),ds))
 def e(i,j):return vals[idx[tuple(sorted((i,j)))]]
 total=A(0)
 for i,j in edges:
  k,l=[a for a in range(4)if a not in (i,j)]
  base=e(i,j);u=e(i,k);v=e(j,k);w=e(i,l);z=e(j,l);opp=e(k,l)
  # Midpoint coordinates reduce interval dependency; equivalent vertex formula retained separately.
  h=u-v;g=w-z
  rr=(u+v)/2-base/4-h*h/(4*base);ss=(w+z)/2-base/4-g*g/(4*base)
  dot=(rr+ss+(h-g)*(h-g)/(4*base)-opp)/2
  assert rr.v>0 and ss.v>0
  t=dot/(rr*ss).sqrt();poly=A(0)
  for c in reversed(P.coeffs()):poly=poly*t+arb(c)
  total+=poly
 return total
stack=[([(Q(1),Q(11,10))]*6,0)];reductions=0;leaves=processed=depth=0;start=time.time()
while stack:
 box,d=stack.pop();processed+=1;depth=max(depth,d);assert d<100 and time.time()-start<240
 z=value(box,True);m=value([((a+b)/2,(a+b)/2)for a,b in box],False).v
 assert z.v.is_finite() and all(v.is_finite()for v in z.d)
 up=m+sum((arb(0,(b-a)/2)*v for (a,b),v in zip(box,z.d)),arb(0))
 if z.v<0 or up<0:leaves+=1;continue
 assert not m>0
 # Uniform derivative signs justify maxima on endpoint faces. Reduce one
 # coordinate at a time, then recompute the complete interval gradient.
 fixed=None
 for j,((lo,hi),g) in enumerate(zip(box,z.d)):
  if lo<hi and (g>0 or g<0):
   end=hi if g>0 else lo;fixed=box.copy();fixed[j]=(end,end);break
 if fixed is not None:
  reductions+=1;stack.append((fixed,d+1));continue
 i=max(range(6),key=lambda i:float(abs(z.d[i]).upper())*float(box[i][1]-box[i][0]));a,b=box[i];mid=(a+b)/2;l=box.copy();r=box.copy();l[i]=(a,mid);r[i]=(mid,b);stack.extend([(l,d+1),(r,d+1)])
out={'pass':True,'dimension':22,'q':'11/10','B_exact':str(P(1)),'strict_upper':'5086.98083741','processed':processed,'leaves':leaves,'depth':depth,'exclusions':0,'monotonic_face_reductions':reductions,'seconds':time.time()-start,'scope':'Auxiliary 3T <= B E only, not a density bound'}
(D/'audit.json').write_text(json.dumps(out,indent=2));print(out)
