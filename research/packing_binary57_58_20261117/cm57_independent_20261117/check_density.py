from pathlib import Path
from fractions import Fraction as F
from math import factorial
import json
P=Path(__file__).resolve().parent
atan=lambda q,N:sum((F((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(N)),F(0))
piL=16*atan(5,200)-4*atan(239,201)
rows=[]
for rec in json.load(open(P/'input/extension_candidates.json')):
 n=rec['dimension'];d=F(rec['center_density_squared']);floor=F({57:'1.53307111021240e-11',58:'5.02419365648403e-12'}[n])
 if n==57:value=d*piL**56*16**29*F(factorial(29),factorial(58))**2
 else:value=d*piL**58/F(factorial(29))**2
 assert floor**2<value
 rows.append(dict(dimension=n,density_floor=str(floor),floor_decimal=float(floor),conditional_on_complete_exclusion=True))
(P/'density_result.json').write_text(json.dumps(rows,indent=2));print(rows)
