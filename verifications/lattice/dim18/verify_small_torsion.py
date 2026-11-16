if not __debug__: raise RuntimeError('Do not run with -O')
"""Independent14D torsion replay by sorted coordinate tuples, not hist recursion."""
from fractions import Fraction as F
from itertools import combinations_with_replacement
import json, math
from pathlib import Path
results={}
for p in [7,9,11]:
 best=-1;arg=None;cases=0
 for v in combinations_with_replacement(range((p+1)//2),18):
  if not any(v):continue
  if p==9 and not any(x%3 for x in v):continue
  cases+=1;costs=[]
  for k in range(1,(p+1)//2):
   residues=[min(k*x%p,p-k*x%p) for x in v]
   numerator=sum(x*x for x in residues)
   if sum(residues)%2:numerator+=p*p-2*p*max(residues)
   costs.append(numerator)
  val=min(costs)
  if val>best:best=val;arg=v
 witness=F(best,2*p*p);assert witness<=F(43,49)
 results[p]={'cases':cases,'physical_squared_witness':str(witness),'coordinates':arg}
 print(p,cases,witness,flush=True)

Path(__file__).with_suffix('.json').write_text(json.dumps(results,indent=2))
