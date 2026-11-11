"""Re-export frozen successful search outputs; does not rerun optimization."""
if not __debug__:raise RuntimeError('Assertions required')
from pathlib import Path
from fractions import Fraction as F
import argparse,json,hashlib,tempfile
from flint import fmpq as Q,fmpq_poly as Poly
D=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--output');args=ap.parse_args()
O=Path(args.output) if args.output else Path(tempfile.mkdtemp(prefix='packing16_exact_export_'))
O.mkdir(parents=True,exist_ok=True)
def write(name,x):(O/name).write_text(json.dumps(x,indent=2)+'\n')
def sh(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bq(x):f=F(float(x));return Q(f.numerator,f.denominator)
for src,dst in [('profile_search/fit.json','fit.json'),('profile_search/inputs/angular_dim16.json','angular.json'),('count_search/factors.json','factors.json')]:
 (O/dst).write_bytes((D/src).read_bytes())
c=json.loads((D/'count_search/reweight_cap0.500300.json').read_text());row=c['history'][-1]
assert row['round']==8 and row['cap']==.5003
write('count_candidate.json',dict(source_scope=c['scope'],factor_sha256=sh(O/'factors.json'),dimension=16,
 cap_exact='5003/10000',rounding='All JSON numeric coefficients are interpreted as their exact printed decimal strings.',record=row))
raw=json.loads((D/'radial_search/conditional_refit.json').read_text())
X=Poly([0,1]);L=[Poly([1]),8-2*X]
for k in range(1,len(raw['radial'])-1):L.append(((2*k+8-2*X)*L[-1]-(k+7)*L[-2])/(k+1))
L=[p/p(0) for p in L]
P=sum((bq(c)*p for c,p in zip(raw['radial'],L)),Poly())
H=sum(((-1)**k*bq(c)*p for k,(c,p) in enumerate(zip(raw['radial'],L))),Poly())
C=bq(raw['C']);lam=bq(raw['lambda_count']);N=H(0);assert N>0
src=json.loads((D/'radial_search/conditional_refit_exact.json').read_text())
assert list(map(str,P.coeffs()))==src['physical'] and list(map(str,H.coeffs()))==src['fourier']
assert str(C)==src['C'] and str(lam)==src['lambda_count_exact']
eps=Q(1,1000000);E=1-X/4
P=P/N+eps*E;H=H/N-eps*E;C/=N;lam/=N
assert H(0)==1-eps
P/=(1-eps);H/=(1-eps);C/=(1-eps);lam/=(1-eps)
write('radial_candidate.json',dict(certified_bound=False,count_conditional=True,source='conditional_refit_exact.json',
 epsilon_exact=str(eps),repair='Normalize source by exact H(0), add epsilon to Laguerre c1, then divide P,H,C,lambda by 1-epsilon.',
 physical=list(map(str,P.coeffs())),fourier=list(map(str,H.coeffs())),C_exact=str(C),lambda_count_exact=str(lam),
 U_exact='7802',control_U_exact='7803',threshold_squared_exact='10000/4997',normalization_exact=str(H(0)),source_normalization_exact=str(N)))
expected={'fit.json':'438c5f0a7d55eca26e6b4c5adbba0c3ec542f0bb18cda15314596456ff2cd4aa',
 'angular.json':'e6f6d85ca4fb160536f3ebc74d94f9ab7ea659793cad94258c8a089447a3f40b',
 'factors.json':'6860eca2c82fdb1b1ebb5d3a4108a216d46c900d6fac7bdeb0d880272e027954',
 'count_candidate.json':'b3225ad4cce819395ce3bd47613f15697af4646fb722ded1916020b5e52ba573',
 'radial_candidate.json':'3c6d8b01a559b52c4418cccca0e93dde7d34118e4cda4a8f398346f91a22be95'}
for name,h in expected.items():assert sh(O/name)==h,('EXPORT HASH MISMATCH',name,sh(O/name),h)
print(json.dumps(dict(status='PASS',scope='Deterministic export from frozen search snapshots, not a new optimization or proof replay.',output=str(O),sha256=expected),indent=2))
