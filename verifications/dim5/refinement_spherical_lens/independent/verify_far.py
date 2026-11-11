"""Independent exact tensor construction and full long-edge lens-domain replay."""
from pathlib import Path
from fractions import Fraction as F
import json,math
import numpy as np
from flint import fmpq as Q,fmpq_poly as P,arb,ctx
if not __debug__:raise RuntimeError('Assertions required')
ctx.prec=384;D=Path(__file__).resolve().parent;s=json.loads((D/'witness.json').read_text());raw=json.loads((D/'witness_psd_exact.json').read_text());budget=next(x for x in json.loads((D/'budget.json').read_text())['dimensions']if x['dimension']==5)
assert (s['dimension'],s['p'],s['K'],F(s['R']))==(5,4,12,F(3,2))
assert [F(str(x))for x in s['x']]==list(map(F,budget['physical_function']['exact_rational_coefficients']))
assert json.loads((D/'budget_audit.json').read_text())['pass']
A=[[[Q(x)for x in row]for row in a]for a in raw['matrices']];bins=[(F(b['lo'])**2,F(b['hi'])**2,F(b['eta']))for b in budget['bins']]
X=P([0,1]);T=[P([1]),2*X-1]
for i in range(1,4):T.append(2*T[1]*T[-1]-T[-2])
# Explicit normalized Gegenbauer coefficients, n=5.
def rising(a,k):return math.prod((a+j for j in range(k)),start=Q(1))
G=[]
for k in range(13):
 p=sum(((-1)**j*rising(Q(3,2),k-j)*(2*X)**(k-2*j)/(math.factorial(j)*math.factorial(k-2*j))for j in range(k//2+1)),P([]));G.append((p/p(1))(-1+Q(16,9)*X))
def bern(p,n):return [sum((p[j]*Q(math.comb(i,j),math.comb(n,j))for j in range(i+1)),Q(0))for i in range(n+1)]
BG=[bern(p,12)for p in G]
def initial(a,b,c,d):
 U=[bern(t(Q(str(a))+Q(str(b-a))*X),4)for t in T];V=[bern(t(Q(str(c))+Q(str(d-c))*X),4)for t in T]
 blocks=[]
 for k in range(13):
  off=int(k==0);blocks.append([[sum((U[i][u]*A[k][i+off][j+off]*V[j][v]for i in range(5)for j in range(5)),Q(0))for v in range(5)]for u in range(5)])
 lo=np.zeros((5,5,13));hi=lo.copy()
 for u,v,w in np.ndindex(lo.shape):
  z=arb(sum((blocks[k][u][v]*BG[k][w]for k in range(13)),Q(0)));lo[u,v,w]=np.nextafter(float(z.lower()),-np.inf);hi[u,v,w]=np.nextafter(float(z.upper()),np.inf)
 return lo,hi

def subdivide(lo,hi,axis):
 lo=np.moveaxis(lo,axis,0);hi=np.moveaxis(hi,axis,0);n=len(lo)-1
 leftlo=[lo[0]];lefthi=[hi[0]];rightlo=[lo[-1]];righthi=[hi[-1]]
 for _ in range(n):
  lo=np.nextafter(np.nextafter(lo[:-1]+lo[1:],-np.inf)/2,-np.inf);hi=np.nextafter(np.nextafter(hi[:-1]+hi[1:],np.inf)/2,np.inf)
  leftlo.append(lo[0]);lefthi.append(hi[0]);rightlo.append(lo[-1]);righthi.append(hi[-1])
 return [np.moveaxis(np.array(z),0,axis)for z in (leftlo,lefthi,rightlo[::-1],righthi[::-1])]
def ranges(box):
 a,b=[1+x/2 for x in box[0]];c,d=[1+x/2 for x in box[1]];t,u=[-1+F(16,9)*x for x in box[2]]
 f=lambda r,s,v:r*r+s*s-2*r*s*v
 maximum=max(f(r,s,t)for r in (a,b)for s in(c,d))
 # Convex quadratic has no interior stationary point in positive quadrant
 # since |u|<1; enumerate minima on each of the four boundary segments.
 pts=[]
 for r in (a,b):pts.append((r,min(d,max(c,u*r))))
 for v in (c,d):pts.append((min(b,max(a,u*v)),v))
 minimum=min(f(r,s,u)for r,s in pts)
 return minimum,maximum,b,d
cut=F(4,25);intervals=[(F(0),cut),(cut,F(1))];stack=[]
for i,(a,b)in enumerate(intervals):
 for c,d in intervals[i:]:
  lo,hi=initial(a,b,c,d);stack.append((lo,hi,[(a,b),(c,d),(F(0),F(1))],0))
passed=excluded=relaxed=processed=depthmax=0
while stack:
 lo,hi,box,depth=stack.pop();processed+=1;depthmax=max(depthmax,depth);assert depth<70 and processed<200000
 assert np.isfinite(lo).all() and np.isfinite(hi).all() and np.all(lo<=hi)
 if hi.max()<0:passed+=1;continue
 wl,wh,r,s=ranges(box)
 if wh<F(9,4):excluded+=1;continue
 allowance=F(0)
 if r<=F(27,25) and s<=F(27,25) and wl>bins[0][0] and wh<bins[-1][1]:
  intersect=[e for a,b,e in bins if not(b<wl or a>wh)];assert intersect;allowance=min(intersect)
 if allowance and hi.max()<np.nextafter(float(arb(Q(str(allowance))).lower()),-np.inf):passed+=1;relaxed+=1;continue
 mid=(lo+hi)/2;scores=[np.max(np.abs(np.diff(mid,axis=i)))for i in range(3)];i=int(np.argmax(scores));ll,lh,rl,rh=subdivide(lo,hi,i);a,b=box[i];m=(a+b)/2;l=box.copy();r=box.copy();l[i]=(a,m);r[i]=(m,b);stack.extend([(rl,rh,r,depth+1),(ll,lh,l,depth+1)])
out={'pass':True,'processed':processed,'passed':passed,'geometric_exclusions':excluded,'relaxed_leaves':relaxed,'maxdepth':depthmax,'scope':'Complete far domain; exact initialization, outward Bernstein subdivision, independently derived distance extrema and conservative closed bin intersections'}
(D/'far_audit.json').write_text(json.dumps(out,indent=2));print(out)
