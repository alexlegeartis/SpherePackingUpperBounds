if not __debug__: raise RuntimeError("Assertions required; do not run with -O")
from fractions import Fraction as F
from itertools import combinations_with_replacement
from pathlib import Path
import json
out={}
for p in [5,7,11,13]:
 best=F(0);cases=0
 for coords in combinations_with_replacement(range(p//2+1),15):
  if not any(coords):continue
  cases+=1;vals=[]
  for k in range(1,p//2+1):
   r=[min(k*x%p,p-k*x%p)for x in coords]
   vals.append(F(sum(x*x for x in r)+(p*p-2*p*max(r))*(sum(r)%2),2*p*p))
  best=max(best,min(vals))
 assert best*F(199,160)<1
 out[p]={'cases':cases,'witness':str(best)};print(p,best,flush=True)
# Independent parity-support average avoids millions of histograms for17,19,23.
for p in [17,19,23]:
 h=p//2;bound=F(0)
 for mask in range(1<<h):
  e=[mask>>j&1 for j in range(h)];m=15 if sum(e)%2 else 14;corr=F(0)
  for k in range(1,h+1):
   rr=[min(k*(j+1)%p,p-k*(j+1)%p)for j in range(h)]
   if sum(a*b for a,b in zip(e,rr))%2:
    corr+=1-F(2*max(r for r,z in zip(rr,e)if z),p)
  bound=max(bound,(F(m*(p+1),12*p)+corr/h)/2)
 assert bound*F(199,160)<1
 out[p]={'parity_supports':1<<h,'witness':str(bound)};print(p,bound,flush=True)
# Even order6: work with coordinates modulo12, reduced by signs to0..6.
# Exact distance computed by nearest-integer choices directly, preserving parity.
best=F(0);cases=0
for coords in combinations_with_replacement(range(7),15):
 if sum(coords)%2:continue
 values=[];exact=True
 for k in range(1,6):
  floors=[k*x//6 for x in coords];res=[k*x%6 for x in coords]
  nearest=[z+(r>3)for z,r in zip(floors,res)]
  if not any(res) and sum(nearest)%2==0:exact=False;break
  distances=[min(r,6-r)for r in res]
  numerator=sum(d*d for d in distances)
  if sum(nearest)%2:numerator+=36-12*max(distances)
  values.append(F(numerator,72))
 if exact:cases+=1;best=max(best,min(values))
assert cases==26720 and best==F(7,9) and best*F(199,160)<1
out['order6']={'cases':cases,'witness':str(best)}
def at(q,N):return sum(((-1)**k*F(1,(2*k+1)*q**(2*k+1))for k in range(N)),F(0))
pi=16*at(5,20)-4*at(239,21)
idxsq=2*F('.0256118267')**2*2027025**2/pi**14*F(151,150)**15
assert idxsq<26**2
assert F(31,8)**2>15 and 1+F(39,32)*F(1,5)==F(199,160)
assert F(1,5)<F(8,14)
allowed=[i for i in range(1,26)if i%6 and all(i%p for p in[5,7,11,13,17,19,23])]
assert max(allowed)==16
out['passed']=True;out['gap']='1/5';out['index_max']=16;out['rho_low']='1024sqrt2'
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2));print('15D THEOREM INDEPENDENTPASS',flush=True)
