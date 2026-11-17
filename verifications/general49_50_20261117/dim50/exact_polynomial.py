import json,math,sys
from fractions import Fraction as F
from flint import fmpq_poly as P,fmpq as Q
from pathlib import Path
sys.set_int_max_str_digits(0);D=Path(__file__).parent
z=json.loads((D/'general_joint_candidate.json').read_text());c=[Q(str(F(str(v))))for v in z['c']];x=P([0,1]);ls=[P([1]),P([25,-2])]
for k in range(2,len(c)):ls.append(((2*k+23-2*x)*ls[-1]-(k+23)*ls[-2])/k)
q=sum((ls[k]*ck*(-1)**k/math.comb(k+24,k)for k,ck in enumerate(c)),P([0]))
(D/'input_roots.json').write_text(json.dumps({'Q':list(map(str,q))}))
