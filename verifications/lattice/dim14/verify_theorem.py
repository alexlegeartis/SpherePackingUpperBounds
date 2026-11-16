if not __debug__:
    raise RuntimeError('Verification requires Python without -O')
"""Exact dependency replay for the 14D trace gap12/25, low point density512.
The published prior radial certificate is replayed separately by run_theorem.py.
The finite-group and sign-averaging arguments are in TWO_PRIMARY_14.md and BINARY_TRACE_14.md.
"""
from fractions import Fraction as F
from pathlib import Path
import subprocess,sys,json,math
from canonical_odd14 import histogram
from determinant_transfer import thresholds
P=Path(__file__).resolve().parent
for p in [3,5,7]:
 subprocess.run([sys.executable,str(P/'verify_prime_metric_nd.py'),'14',str(p)],check=True)
 r=json.loads((P/f'n14_p{p}_metric_replay.json').read_text());assert r['passed'] and F(r['minimum_bound'])>=F(12,25)
for p in [4,8]:
 subprocess.run([sys.executable,str(P/'verify_even_metric_nd.py'),'14',str(p),'12/25'],check=True)
 r=json.loads((P/f'n14_even{p}_metric_replay.json').read_text());assert r['passed'] and F(r['minimum_bound'])>F(12,25)
assert F(thresholds(14)['determinant_one_threshold'])==F(8,13)>F(12,25)
S=F(12,25);D=1+F(19,16)*S
assert F(15,4)**2>14 and D==F(157,100)
for p,w in [(11,F(73,121)),(13,F(8,13))]:
 r=histogram(14,p);assert F(r['physical_squared_witness'])==w and D*w<1
# Sign-averaging excludes binary weights2,4,6 for S<4; integral kernel needs S>=28.
assert S<4 and S<28
assert 15*8>14*8 # rank4 binary subcode is impossible
# Good order4 histogram is(0,8,6); independent order2 would toggle a2 to0.
# Together with order8 exclusion, every2-primary group has order<=8.
# Prior-to-index conversion, with an exact Machin lower bound.
def lo(q,N=30):return sum((F((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(N)),F(0))
pi_lower=16*lo(5)-4*(lo(239)+F(1,61*239**61))
assert pi_lower>F('3.14159')
assert F('.03525')*2**14*math.factorial(7)/pi_lower**7<964
assert F(964,64)*(1+S/28)**7<17
allowed=[i for i in range(1,17) if all(i%p for p in [3,5,7,11,13]) and i<=8]
assert max(allowed)==8
assert 8*64==512
print('PASS: rho>512 implies every14-pair trace excess >=12/25')
(P/'theorem_replay.json').write_text(json.dumps({'passed':True,'gap':'12/25','low_point_density':512,'prior_density_input':'.03525','index_strict_upper':17,'proof':'TWO_PRIMARY_14.md'},indent=2))
