from pathlib import Path
from fractions import Fraction as F
import json,math
from flint import fmpq as Q,fmpq_poly as P
if not __debug__:raise RuntimeError('Assertions required')
d=Path(__file__).resolve().parent;g=json.loads((d/'kernel_guard5.json').read_text());cs=list(map(Q,g['gegenbauer_coefficients']));assert g['dimension']==5 and cs[0]==Q(1,2) and min(cs)>=0
X=P([0,1]);lam=Q(3,2)
def rising(a,k):return math.prod((a+j for j in range(k)),start=Q(1))
basis=[]
for k in range(len(cs)):
 p=sum(((-1)**j*rising(lam,k-j)*(2*X)**(k-2*j)/(math.factorial(j)*math.factorial(k-2*j))for j in range(k//2+1)),P([]));basis.append(p/p(1))
p=sum((a*b for a,b in zip(cs,basis)),P([]))+Q(1,2);seq=[p,p.derivative()]
while seq[-1].degree()>0:
 z=-(seq[-2]%seq[-1]);assert z;seq.append(z)
def V(t):
 a=[int(p(t)>0)-int(p(t)<0)for p in seq if p(t)];return sum(x!=y for x,y in zip(a,a[1:]))
assert p(-1)<0 and p(Q(7,9))<0 and V(Q(-1))==V(Q(7,9))
D=sum(cs);a=D*D/(4*cs[0]);b=-D/2;assert a*cs[0]==b*b and 2*b+D==0
assert a==Q(g['kernel_origin']) and b==Q(g['kernel_cross']) and D==Q(g['kernel_diagonal'])
print('PASS independent exact guard PSD, zero pair correction and strict angular sign')
