import sys,json,math
if not __debug__: raise RuntimeError("Assertions must be enabled")
sys.set_int_max_str_digits(0);sys.path.insert(0,'solver_packages')
from fractions import Fraction as F
from flint import fmpq_poly,fmpq
from pathlib import Path
p=Path(__file__).resolve().parent;ff=lambda x:fmpq(x.numerator,x.denominator)
def rise(x,k):return math.prod(x+i for i in range(k))
rows=[]
for name in ['cyclic_angular40_exact.json']:
 d=json.loads((p/name).read_text());cs=list(map(F,d['gegenbauer_coefficients']));ps=list(map(F,d['polynomial_coefficients']));assert cs[0]==1 and min(cs)>=0 and not any(cs[20:]);cs=cs[:20]
 got=[F(0)]*20
 for l,b in enumerate(cs):
  norm=F(rise(F(d['dimension']-2),l),math.factorial(l))
  for j in range(l//2+1):got[l-2*j]+=b*(-1)**j*rise(F(d['dimension']-2,2),l-j)*2**(l-2*j)/(math.factorial(j)*math.factorial(l-2*j)*norm)
 assert got==ps and sum(cs)==F(d['B_exact'])
 poly=fmpq_poly(list(map(ff,ps)));seq=[poly,poly.derivative()]
 while seq[-1].degree()>0:
  rem=-(seq[-2]%seq[-1]);assert rem
  num=rem.numer();seq.append(fmpq_poly(num/abs(num.content())))
 def var(t):
  vals=[g(ff(t)) for g in seq];sg=[int(v>0)-int(v<0) for v in vals if v];return sum(a!=b for a,b in zip(sg,sg[1:]))
 cut=F(d['cutoff']);assert poly(-1)<0 and poly(ff(cut))<0
 v0=var(F(-1));v1=var(cut);assert v0==v1;assert var(cut)-var(F(1))==1 and poly(1)>0
 cap=poly+fmpq(1,100000);seq=[cap,cap.derivative()]
 while seq[-1].degree()>0:
  rem=-(seq[-2]%seq[-1]);num=rem.numer();seq.append(fmpq_poly(num/abs(num.content())))
 assert var(F(-1))==var(F(1,2)) and cap(-1)<0 and cap(fmpq(1,2))<0
 rows.append({'file':name,'pass':True,'sturm_variations':[v0,v1],'B':str(sum(cs))})

(p/'independent_cyclic40_angular.json').write_text(json.dumps(rows,indent=2));print(rows)
