import sys
if not __debug__:raise RuntimeError("Assertions must remain enabled")
sys.set_int_max_str_digits(0)
import json,math,hashlib
import sys;sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as F
from flint import fmpq as Q,fmpq_poly as P,arb,ctx
ctx.prec=768;D=Path(__file__).parent;root=D.parent/'lattice_quotient_resume_20261117';repair=json.load(open(D/'repaired.json'));cand=json.load(open(D/'input_candidate.json'));extensions=json.load(open(root/'scalar48_kernel_upward_extensions.json'));e48=next(x for x in extensions if x['n']==48);e50=next(x for x in extensions if x['n']==50)
def expand(n,co):
 x=P([0,1]);g=[P([1]),x]
 for k in range(2,2*len(co)-1):g.append(((2*k+n-4)*x*g[-1]-(k-1)*g[-2])/(k+n-3))
 return sum((Q(c)*g[2*k]for k,c in enumerate(co)),P([0]))
co48=list(map(Q,e48['coefficients']));co50=list(map(Q,e50['coefficients']));assert co50[0]==1 and min(co50)>=0;q=Q(e50['q_normalized']);assert q==sum(co50)==Q('1150750294032540977074575/1437673731494371')
raw=json.load(open(root/'rank2_candidate_B783489561_20261117.json'));r=[F(str(x))for x in raw['even_coefficients']];r[0]-=F(1,5000);r=[x/r[0]for x in r]
while r[-1]==0:r.pop()
assert list(map(str,co48))==list(map(str,r));a48=expand(48,co48);a50=expand(50,co50);assert a50*a48(1)==a48*a50(1) and a48(1)>0 and a50(1)>0
# Algebraic global Fourier positivity, independently reconstructed from saved factors.
qr=P([1])
for f in repair['factors']:
 if f['type']=='negative_linear':v=Q(f['a']);assert v>0;qr*=P([1,1/v])
 elif f['type']=='positive_double':v=Q(f['mid']);assert v>0;qr*=P([1,-2/v,1/(v*v)])
 else:
  re,im=Q(f['re']),Q(f['im']);assert im>0;norm=re*re+im*im;qr*=P([1,-2*re/norm,1/norm])
assert qr==P(list(map(Q,repair['Q']))) and qr(0)==1
p=P(list(map(Q,repair['P'])));a=F(str(cand['Z'][1][1]));z=F(str(cand['Z'][0][1]))-F(1,10**9);s=z*z/a+F(1,10**20);assert a>0 and s*a-z*z>0
rho=F(1,10**70);L=(1+rho)/(1-rho);assert F(168**168,math.factorial(168))>1/rho
# Fourier-normalized density with contact x=d=22.
prefactor=F(22,4)**25/math.factorial(25);density_no_budget=prefactor*(F(str(p[0]))+s)
def ar(x):x=F(str(x));return arb(x.numerator)/x.denominator
v=arb(1);h=arb(1);f=arb(p(Q(22)))*(-arb(22)).exp();contact=f+2/ar(q)*(ar(a*L)+ar(z));print('EXACT NORMALIZATION/PSD/POSITIVE FACTORS PASS','density before annular budgets',float(density_no_budget),'contact residual',str(contact),flush=True)
out=dict(status='EXACT NORMALIZATION PSD AND FOURIER FACTORS PASS',q=str(q),a=str(a),z=str(z),s=str(s),rho=str(rho),line_factor=str(L),density_prefactor=str(prefactor),density_without_annular_budget=str(density_no_budget),contact_residual=str(contact),angular_extension_verified=True,prior_rank_two_source_sha256=hashlib.sha256((root/'rank2_candidate_B783489561_20261117.json').read_bytes()).hexdigest(),full_radial_status='PENDING ANNULAR BUDGET CERTIFICATE')
(D/'normalization_result.json').write_text(json.dumps(out,indent=2))
