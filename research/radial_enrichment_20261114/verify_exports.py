"""Reconstruct the three frozen rational exports and check family extensions."""
import json,sys,hashlib
from fractions import Fraction
from pathlib import Path
from flint import fmpq as Q, fmpq_poly as P, arb_poly, arb, ctx
from model import inputs,ROOT,angular_path,pi_bounds

if not __debug__:raise RuntimeError('Assertions must remain enabled')
sys.set_int_max_str_digits(0);ctx.prec=512
HERE=Path(__file__).resolve().parent
rows=[]
pl,pu=pi_bounds()
for n,d,previous_degree in [(45,47,31),(46,63,47),(47,47,47)]:
    package=ROOT/f'verifications/dim{n}/refinement_radial20261114'
    r=json.loads((package/'certificate.json').read_text())
    cand=HERE/f'n{n}_d{d}_c40/candidate.json'
    c=json.loads(cand.read_text());assert c['dimension']==n and c['degree']==d
    cs=[Q(str(Fraction(float(v)).limit_denominator(10**14))) for v in c['x']]
    source,h,ls,B=inputs(n,d)
    assert len(cs)==d+3 and cs[0]>=0 and cs[-1]>0
    correction=sum((cs[k+1]*ls[k] for k in range(d+1)),P([]))
    pc=cs[0]*source+correction
    hc=cs[0]*h+sum(((-1)**k*cs[k+1]*ls[k] for k in range(d+1)),P([]))
    assert pc==P(list(map(Q,r['physical_polynomial'])))
    assert hc==P(list(map(Q,r['fourier_polynomial'])))
    assert Q(r['M_exact'])==cs[-1]/B*Q(1000001,1000000)
    assert Q(r['B_exact'])==B and correction.degree()==d and source.degree()==111
    assert pc.degree()==hc.degree()==111
    added=[k for k in range(previous_degree+1,d+1) if cs[k+1]]
    if n in [45,46]:assert added and correction.degree()>previous_degree
    ap=arb_poly([arb(v) for v in pc.coeffs()]);escape=None
    # Exact rational points, followed by exact evaluation; sampling chooses a witness only.
    for j in range(2200,4000):
        x=Q(j,100)
        if x<=7*pu:continue
        if ap(arb(x))>0:
            exact=pc(x);assert exact>0 and 7*pu<x<40
            escape={'point_X_exact':str(x),'physical_polynomial_value_exact':str(exact),
                    'positive_exact':True,'beyond_old_7pi_exact':True}
            break
    assert escape is not None
    (package/'tail_escape.json').write_text(json.dumps(escape,indent=2)+'\n')
    rows.append({'dimension':n,'degree':d,'source_t_nonnegative':True,
                 'physical_export_exact':True,'Fourier_export_exact':True,'M_export_exact':True,
                 'correction_degree_exact':d,'previous_correction_degree':previous_degree,
                 'nonzero_new_Laguerre_degrees':added,'old_7pi_tail_template_violated':True,
                 'candidate_sha256':hashlib.sha256(cand.read_bytes()).hexdigest(),
                 'certificate_sha256':hashlib.sha256((package/'certificate.json').read_bytes()).hexdigest(),
                 'angular_sha256':hashlib.sha256(angular_path(n).read_bytes()).hexdigest()})
(HERE/'EXPORT_CHECKS.json').write_text(json.dumps(rows,indent=2)+'\n')
print('PASS exact exports and explicit positive old-tail witnesses for 45D,46D,47D')
