"""Exact dimension-transfer constants for same-mod-2-class frame determinants.
No optimizer or floating-point computation is used for theorem verification.
"""
from fractions import Fraction as F
import json

def thresholds(n):
 assert isinstance(n,int) and n>=3
 star=F(8,n-1)
 negative=F(4*(n-2)**2,3*n*n-8*n+6)
 safe=min(star,negative)
 # At safe, all vertex matrices are positive definite: star via
 # S(sqrt(n)-1)<4, verified without irrational arithmetic.
 assert (1+F(4,safe))**2>n
 z=safe/(n-2)
 assert (n-1)*z<2
 assert F(n-2)*z/2-F((n-1)**2+1)*z*z/(8*(1-F(n-1)*z/2))>=0
 assert 1+safe/2-F(n-1,16)*safe*safe>=1
 return dict(n=n,star_threshold=str(star),negative_vertex_threshold=str(negative),determinant_one_threshold=str(safe),covolume_baseline_squared=str(F(2)**(2-n)))

if __name__=='__main__':
 rows=[thresholds(n) for n in range(3,65)]
 assert all(F(r['determinant_one_threshold'])==F(r['star_threshold']) for r in rows if r['n']>=9)
 out={'theorem':'min lattice length1; n distinct antipodal pairs in one mod2 class; includes one short pair','rows':rows,'passed':True}
 from pathlib import Path
 Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2))
 for r in rows:
  if r['n'] in [3,4,5,6,7,8,9,12,16,20,24,32,48,64]: print(r)
