import json,math,sys
from pathlib import Path
from fractions import Fraction as F
sys.set_int_max_str_digits(0)
if not __debug__:raise RuntimeError('Assertions required')
D=Path(__file__).parent
z=json.load(open(D/'root_repair44/annular_final_certificate.json'))
value=(F(z['P'][0])+F(z['s'])+F(z['budget']))*F(22,4)**22/math.factorial(22)
safe=F('0.000001776548487')
assert 0<value<safe
(D/'density_result.json').write_text(json.dumps(dict(status='PRODUCER COMPLETE; INDEPENDENT AUDIT PENDING',safe_upper=str(safe),density_exact=str(value)),indent=2))
print(float(value),float(safe))
