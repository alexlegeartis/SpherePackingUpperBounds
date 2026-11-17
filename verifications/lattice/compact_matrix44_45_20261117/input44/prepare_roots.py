import json,sys
from pathlib import Path
from fractions import Fraction as F
from flint import fmpq_poly as P,fmpq as Q
sys.set_int_max_str_digits(0);D=Path(__file__).parent/'root_repair44';z=json.load(open(D/'input_candidate.json'));assert z['n']==44;c=list(map(lambda x:F(str(x)),z['c']));alpha=F(21);x=P([0,1]);L=[P([1]),P([22,-2])]
for k in range(2,len(c)):L.append(((Q(str(2*k+alpha-1))-2*x)*L[-1]-Q(str(k+alpha-1))*L[-2])/k)
q=P([0]);norm=F(1)
for k,v in enumerate(c):
 if k:norm*=F(k+alpha,k)
 q+=L[k]*Q(str(v*(-1)**k/norm))
(D/'input_roots.json').write_text(json.dumps(dict(Q=list(map(str,q)))))
