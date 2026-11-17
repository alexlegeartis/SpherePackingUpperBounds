import math
from flint import arb
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
 def __pow__(self,k):
  assert isinstance(k,int) and k>=0
  if k==0:return Jet(1)
  if k==1:return self
  v=self.v**k;d=k*self.v**(k-1);dd=k*(k-1)*self.v**(k-2);return Jet(v,[d*x for x in self.g],[[d*self.h[i][j]+dd*self.g[i]*self.g[j]for j in range(3)]for i in range(3)])
 def log(self):return Jet(self.v.log(),[x/self.v for x in self.g],[[self.h[i][j]/self.v-self.g[i]*self.g[j]/self.v**2 for j in range(3)]for i in range(3)])
 def exp(self):
  v=self.v.exp();return Jet(v,[v*x for x in self.g],[[v*(self.h[i][j]+self.g[i]*self.g[j])for j in range(3)]for i in range(3)])
def pe2(cs,z):
 mid=z.v.mid();delta=z.v-mid;n=len(cs);tt=list(cs)
 for i in range(n-2,-1,-1):
  for j in range(i,n-1):tt[j]=tt[j]+mid*tt[j+1]
 value=tt[-1];der=arb(0);sec=arb(0)
 for c in tt[-2::-1]:sec=sec*delta+2*der;der=der*delta+value;value=value*delta+c
 return Jet(value,[der*x for x in z.g],[[sec*z.g[i]*z.g[j]+der*z.h[i][j]for j in range(3)]for i in range(3)])
