"""Read-only numerical/structural QA and report generation; no additional fit/LP."""
import hashlib
import json
from fractions import Fraction

import numpy as np

from log_hinge import D, V_KNOTS, LOG_KNOTS, LogHingeExperiment, write

e=LogHingeExperiment(16)
fit=json.loads((D/'fit.json').read_text())
radial=json.loads((D/'radial_result.json').read_text())
p=np.array(fit['parameters'])
controls=json.loads((D/'frozen_control/dim16_dual_refit_radial_comparison.json').read_text())['outcomes']
control=next(x for x in controls if x['index']=='dual_refit')
quadratic=next(x for x in controls if x['index']=='quadratic')
old=np.r_[json.loads((D/'frozen_control/dim16_dual_refit.json').read_text())['candidate']['parameters'],np.zeros(7)]

# The first seven feature arrays are byte-for-byte the old v-hinge arrays.
v=np.unique(np.r_[1.,V_KNOTS,LOG_KNOTS,np.linspace(1,32,2001)])
oldbasis=np.c_[v-1,np.maximum(v[:,None]-V_KNOTS,0)]
assert np.array_equal(e.basis(v)[:,:7],oldbasis)
exact_parameters=list(map(lambda x:Fraction(str(x)),fit['parameters']))
assert min(exact_parameters)>=0
final_slope=sum(exact_parameters[:7])
assert final_slope>=7

# Analytic profile-parameter Jacobian agrees with finite differences.
rows=e.rows(np.array([-.8,.181,.35,.2,-.15]),np.array([1.1,1.,1.2,1.5,2.3]))
_,J=e.expression(rows,p,jac=True)
error=0.
for k in range(14):
    step=np.zeros(14); step[k]=1e-5
    fd=(e.expression(rows,p+step)-e.expression(rows,p-step))/(2e-5)
    error=max(error,float(np.max(np.abs(fd-J[:,k]))))
assert error<1e-8

prices=np.load(D/'candidate_pair_dual_pricing.npz')
x=prices['pair'];v=x/(2*np.pi)
h=np.exp(-e.basis(v)@p)
weights=(-prices['scaled_marginals']/prices['fixed_row_scale'])*2*radial['C']/e.B*np.exp(x)*h*(1-h)
active=np.where(weights>weights.max()*1e-7)[0]
def obstacle(v):
    q0=e.q(1-v/2);q2=e.q(np.sqrt(v)/2)
    assert np.all(q0<0) and np.all(q2>0)
    return np.log(2*q2/-q0)
def profile_row(v,weight=None):
    vv=np.array([v]); lower=float(obstacle(vv)[0])
    newphi=float((e.basis(vv)@p)[0]); oldphi=float((e.basis(vv)@old)[0])
    out=dict(v=float(v),isosceles_lower_obstacle=lower,phi=newphi,
        slack=newphi-lower,old_control_phi=oldphi,old_control_slack=oldphi-lower)
    if weight is not None:out['current_pair_tangent_weight']=float(weight)
    return out

summary=dict(certified=False,status='bounded experiment completed; no promotion',
    source_sha256=e.sha256,source_Q_unchanged=True,old_seven_feature_arrays_exactly_preserved=True,
    analytic_jacobian_max_absolute_error=error,nonnegative_parameters=True,
    final_linear_slope_exact_decimal=str(final_slope),
    Gaussian_tail_floor=7,
    sampled_expression=radial['density_expression'],
    old_cone_matched_expression=control['density_expression'],
    quadratic_matched_expression=quadratic['density_expression'],
    improvement_over_old_cone=control['density_expression']-radial['density_expression'],
    comparator=radial['comparator'],gap_to_comparator=radial['gap_to_comparator'],
    candidate_dual_rows=[profile_row(v[i],weights[i]) for i in active],
    old_primary_dual_rows=[profile_row(z) for z in [1.61532400564816,1.7298028904199105]],
    full_proof_attempted=False,compact_cyclic_gate='NUMERICAL ONLY',
    Fourier_sign_gate='NOT RUN',physical_tail_sign_gate='NOT RUN',
    annular_radial_gate='NOT RUN',denied_actions=[],blockers=[])
write('summary.json',summary)
manifest={str(path.relative_to(D)):hashlib.sha256(path.read_bytes()).hexdigest()
          for path in sorted(D.rglob('*')) if path.is_file() and '__pycache__' not in str(path)
          and path.name not in ['manifest.json','summary.log']}
write('manifest.json',manifest)
print(json.dumps(summary,indent=2))
