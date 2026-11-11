"""Replay embedded degree31 LOWER certificates and compare exact thresholds."""
from fractions import Fraction as F
from pathlib import Path
import json
import subprocess
import sys

if not __debug__:
    raise RuntimeError('Do not use Python -O')
D=Path(__file__).resolve().parent
C=D/'lower_comparison'
radial=json.loads((D/'radial_result.json').read_text())
n=radial['dimension']
assert radial['pass'] and radial['outside_degree31_span_exact']
upper=F(radial['strict_density_upper'])
angular={}
for tag,filename in [('joint',f'angular{n}_lower_certificate.json'),('allrate',f'angular{n}_caponly_certificate.json')]:
    path=C/'inputs'/filename
    if tag=='allrate':
        data=json.loads(path.read_text())
        assert all(row['point'][0]=='cap' for row in data['dual_rows'])
    output=C/f'angular_{tag}_replayed.json'
    subprocess.run([sys.executable,str(C/'verify_angular_certificate.py'),str(path),'--out',str(output)],check=True)
    angular[tag]=json.loads(output.read_text())
    assert angular[tag]['pass'] and angular[tag]['all_harmonic_degrees_certified']

floors={}
saved_B=F(json.loads((D/f'inputs/angular{n}.json').read_text())['B_exact'])
for tag in ['fixed','joint','allrate']:
    selected=C/f'dual{n}_{tag}_selected.json'
    subprocess.run([sys.executable,str(C/'verify_dual.py'),str(selected)],check=True)
    record=json.loads((C/f'dual{n}_{tag}_verified.json').read_text())
    assert record['pass'] and record['all_pair_rows_at_contact']
    if tag=='fixed':assert F(record['B_exact'])==saved_B
    else:assert F(record['B_exact'])==F(angular[tag]['lower_bound_B_exact'])
    lower=F(record['safe_density_lower_exact'])
    floors[tag]={'lower_exact':str(lower),'lower_decimal':record['safe_density_lower_decimal'],
                 'new_upper_below_floor':upper<lower,
                 'floor_minus_new_upper_exact':str(lower-upper)}

assert floors['fixed']['new_upper_below_floor']
if n==46:
    assert floors['allrate']['new_upper_below_floor'] and floors['joint']['new_upper_below_floor']
elif n==47:
    assert not floors['joint']['new_upper_below_floor'] and not floors['allrate']['new_upper_below_floor']
else:raise RuntimeError('Unexpected dimension')

# Exact universal amplitude identity, not a sampled numerical inequality:
# A^2-D(2A-D)=(A-D)^2. Combined with CD>=B*A^2 and D>0,
# 2A-D<=A^2/D<=C/B. For D=0, PSD forces A=0.
out={'pass':True,'dimension':n,'new_degree47_strict_upper':radial['strict_density_upper'],
     'all_embedded_lower_proofs_replayed':True,'degree31_floors':floors,
     'outside_degree31_span_exact':True,
     'free_scalar_amplitude_extension':'For B>0,C,D>=0, CD>=B*A^2 and contact f<=2A-D imply f<=C/B; D=0 is handled separately.',
     'scope':'Degree31 lower floors retain the normalized source span and physical negativity beyond7pi; allrate permits every positive scalar Gaussian rate and all finite angular degrees. Harmonic-dependent radial matrices are excluded.'}
(D/'saturation_comparison_result.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
