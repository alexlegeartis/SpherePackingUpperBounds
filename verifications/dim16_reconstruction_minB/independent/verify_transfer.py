"""Audit exact algebraic transfer; baseline radial/count requires its separate replay."""
if not __debug__: raise RuntimeError('Assertions required')
from fractions import Fraction as F
from pathlib import Path
import json,math,hashlib
import sys; sys.set_int_max_str_digits(0)
D=Path(__file__).resolve().parent
old=D.parents[1]/'dim16'
new=D/'geometry'
rd=lambda p:json.loads(p.read_text())
a=rd(old/'geometry/angular.json');b=rd(new/'angular.json');r=rd(old/'radial/candidate_exact.json')
assert rd(new/'VERIFIED.json')['status']=='PASS'
for p,h in rd(new/'VERIFIED.json')['input_hashes'].items(): assert hashlib.sha256((new/p).read_bytes()).hexdigest()==h
assert (old/'geometry/fit.json').read_bytes()==(new/'fit.json').read_bytes()
assert (old/'radial/inputs/fit.json').read_bytes()==(new/'fit.json').read_bytes()
Bo,Bn=F(a['B_exact']),F(b['B_exact']);Co=F(r['C_exact']);Cn=Co*Bn/Bo
assert 0<Bn<Bo and 0<Cn<Co and Cn/Bn==Co/Bo
assert F(r['normalization_exact'])==1 and F(r['U_exact'])==7802
assert F(r['threshold_squared_exact'])==F(10000,4997)
# No pair function or count allowance changes at all.
def arctan(q,N=160):
 s=sum((F((-1)**k,(2*k+1)*q**(2*k+1))for k in range(N)),F(0))
 t=s+F((-1)**N,(2*N+1)*q**(2*N+1));return min(s,t),max(s,t)
a,b=arctan(5);c,d=arctan(239);pl,pu=16*a-4*d,16*b-4*c
objective=F(r['physical'][0])+Cn+7802*F(r['lambda_count_exact'])
lo=objective*pl**8/(256*math.factorial(8));hi=objective*pu**8/(256*math.factorial(8))
assert 0<lo<hi<F('0.024920883')
out={'status':'PASS','scope':'Exact transfer plus separately audited geometry. Requires baseline original radial/count replays at commit043b119.','bound':'0.024920883','C_exact':str(Cn),'density_upper_rational':str(hi),'density_display':float(hi),'same_pair_function':True,'new_geometry_audit':'geometry/VERIFIED.json'}
(D/'TRANSFER_AUDIT.json').write_text(json.dumps(out,indent=2)+'\n');print(out['status'],out['density_display'])
