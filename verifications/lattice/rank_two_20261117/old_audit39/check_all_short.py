if not __debug__:raise RuntimeError('Assertions required')
from fractions import Fraction as R
from itertools import combinations
import json
from pathlib import Path
results=[]
for dim in (39,):
 tested=skipped=0
 for tick in range(1,97):
  delta=R(tick,384);slope=R(21,20)*(dim-1)*delta/(4*(1+delta))
  c=delta/(2+delta);whole=dim*(1-c*c)/(1-dim*c*c)
  for count in range(dim+1,int(whole)+1):
   # Avoid floored A formula: relative cap >=count-1 is equivalent to this square comparison.
   threshold=R(count-dim,(dim-1)*(count-2))
   def permitted(z):
    eta=(delta-z)**2/(4*(2-z)*(2+delta));cc=(c+eta)/(1-eta)
    return cc*cc>=threshold
   if not permitted(R(0)):
    skipped+=1;continue
   zbound=delta
   if not permitted(delta):
    lower=R(0);upper=delta
    for _ in range(16):
     middle=(lower+upper)/2
     if permitted(middle):lower=middle
     else:upper=middle
    zbound=upper
   limit=R(count-dim)/slope
   lines=[(R(1),R(0),R(0)),(R(1),R(0),zbound),(R(2-count),R(1),R(0)),(R(1),R(1),(count-1)*delta),(R(0),R(1),limit)]
   vertices=set()
   for (a,b,e),(f,g,h) in combinations(lines,2):
    determinant=a*g-f*b
    if determinant:
     z=(e*g-h*b)/determinant;S=(a*h-f*e)/determinant
     if 0<=z<=zbound and (count-2)*z<=S<=(count-1)*delta-z and S<=limit:vertices.add((z,S))
   for z,S in vertices:
    second=(delta+z)*S+2*z*z-(count-2)*delta*z
    value=dim*(4*count+4*S+R(count,2)*second+S*S/2)-(2*count+S)**2
    assert value<=0,(dim,tick,count,z,S,value)
    assert S==limit or value<0,(dim,tick,count,'strict')
   for a,b,e in lines:
    dz,dS=b,-a
    assert dim*count*dz*dz+R(dim*count,2)*dz*dS+(R(dim,2)-1)*dS*dS>=0
   tested+=1
  assert 1+slope<dim # singleton, z<=1
 results.append(dict(n=dim,cutoffs=96,checked_cases=tested,impossible_by_projected_cap=skipped,pass_exact=True))
print(json.dumps(results,indent=2))
Path(__file__).with_suffix('.json').write_text(json.dumps(results,indent=2))
