"""Portable full proof of the combined 16D density upper bound."""
if not __debug__:
    raise RuntimeError('Assertions are required; do not use -O or -OO.')
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,os,subprocess,sys,time
from flint import arb,fmpq,ctx
ctx.prec=256
D=Path(__file__).resolve().parent
C=D/'count/density_kissing16_independent_cyclic_20261111'
S=D/'count/density_kissing16_reweight_frozen_20261111'
G=D/'geometry';R=D/'radial'
EXPECTED={
 'geometry/fit.json':'438c5f0a7d55eca26e6b4c5adbba0c3ec542f0bb18cda15314596456ff2cd4aa',
 'geometry/angular.json':'e6f6d85ca4fb160536f3ebc74d94f9ab7ea659793cad94258c8a089447a3f40b',
 'radial/candidate_exact.json':'3c6d8b01a559b52c4418cccca0e93dde7d34118e4cda4a8f398346f91a22be95',
 'count/density_kissing16_reweight_frozen_20261111/candidate.json':'b3225ad4cce819395ce3bd47613f15697af4646fb722ded1916020b5e52ba573',
 'count/density_kissing16_reweight_frozen_20261111/factors.json':'6860eca2c82fdb1b1ebb5d3a4108a216d46c900d6fac7bdeb0d880272e027954'}
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def a(q):q=F(q);return arb(fmpq(q.numerator,q.denominator))
for p,h in EXPECTED.items():assert sha(D/p)==h,('HASH',p)
assert (G/'fit.json').read_bytes()==(R/'inputs/fit.json').read_bytes()
assert (G/'angular.json').read_bytes()==(R/'inputs/angular_dim16.json').read_bytes()
env=os.environ.copy()
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:env[key]='1'
env.pop('PYTHONOPTIMIZE',None)
stages=[];start=time.monotonic();logs=D/'logs';logs.mkdir(exist_ok=True)
for name,path in [('spherical_count',C/'audit_count.py'),('angular_geometry',G/'verify_geometry.py'),('radial_pair',R/'verify.py')]:
 t=time.monotonic();run=subprocess.run([sys.executable,str(path)],cwd=D,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
 (logs/(name+'.log')).write_text(run.stdout)
 assert run.returncode==0,(name,run.returncode)
 stages.append(dict(stage=name,status='PASS',seconds=time.monotonic()-t));print('PASS',name,flush=True)
count=read(C/'CERTIFIED_COUNT.json');count_run=read(C/'AUDIT_RESULT.json')
geo=read(G/'RESULT_GEOMETRY.json');rad=read(R/'RESULT.json');candidate=read(R/'candidate_exact.json')
assert count_run['status']=='PASS' and geo['status']=='PASS' and rad['status']=='PASS'
U=F(count['real_upper_bound_exact']);integer_U=U.numerator//U.denominator
assert integer_U==7802==count['spherical_count_integer_upper_bound']==rad['count_U']
assert F(candidate['U_exact'])==integer_U
cap=F(count['cap_exact']);tau=F(candidate['threshold_squared_exact']);s=F(2)
assert cap==F(5003,10000) and tau==F(10000,4997) and s<tau<2*s
assert 1-s/(2*tau)==cap
P=list(map(F,candidate['physical']));H=list(map(F,candidate['fourier']))
C0=F(candidate['C_exact']);lam=F(candidate['lambda_count_exact'])
assert C0>=0 and lam>=0 and H[0]>0
value=arb.pi()**8/40320*a((s/4)**8)*a((P[0]+C0+integer_U*lam)/H[0])
target=F('0.024920988');old=F('0.02492121');lower=arb.pi()**8/645120
assert value<a(target)<a(old)
out=dict(status='PASS',dimension=16,strict_density_upper_exact=str(target),strict_density_upper_decimal='0.024920988',
 density_enclosure=str(value),complete_density_claim=True,
 spherical_count=integer_U,spherical_cap_exact=str(cap),minimum_squared_distance=str(s),
 half_open_shell_threshold_squared_exact=str(tau),
 same_profile_and_angular_inputs=True,all_component_gates_replayed=True,
 comparator_exact=str(old),achieved_lower_formula='pi^8/645120',achieved_lower_enclosure=str(lower),
 upper_reduction_percent=str(100*(a(old)-a(target))/a(old)),
 gap_closure_percent=str(100*(a(old)-a(target))/(a(old)-lower)),
 stages=stages,seconds=time.monotonic()-start,input_sha256=EXPECTED,
 scope='Complete combined density certificate. Separate independent radial/geometry audit is distributed alongside this primary package.')
(D/'RESULT.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2),flush=True)
