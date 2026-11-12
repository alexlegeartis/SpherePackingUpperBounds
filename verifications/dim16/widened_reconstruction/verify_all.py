"""Portable complete geometry, radial and spherical count replay."""
if not __debug__:raise RuntimeError('Assertions required')
from pathlib import Path
from fractions import Fraction as F
import sys,subprocess,json,hashlib
D=Path(__file__).resolve().parent
components_only='--components-only' in sys.argv
read=lambda p:json.loads((D/p).read_text())
for p,h in read('INPUT_HASHES.json').items():assert hashlib.sha256((D/p).read_bytes()).hexdigest()==h,('HASH',p)
for f,r in [('fit.json','fit.json'),('angular.json','angular_dim16.json')]:
 data=(D/'geometry'/f).read_bytes()
 for p in [D/'geometry/independent'/f,D/'radial/inputs'/r,D/'radial/independent/inputs'/r]:assert data==p.read_bytes(),('LINKAGE',str(p))
for f in ['candidate_exact.json','conditional_refit.json','conditional_refit_exact.json']:
 assert (D/'radial'/f).read_bytes()==(D/'radial/independent'/f).read_bytes(),('RADIAL IDENTITY',f)
for f in ['geometry/verify_geometry.py','geometry/independent/verify.py','radial/verify.py','radial/independent/verify.py']:subprocess.run([sys.executable,str(D/f)],check=True)
assert read('geometry/RESULT_GEOMETRY.json')['status']=='PASS'
assert read('geometry/independent/VERIFIED.json')['status']=='PASS'
r=read('radial/candidate_exact.json');ri=read('radial/independent/VERIFIED.json');rr=read('radial/RESULT.json')
assert ri['status']==rr['status']=='PASS' and F(r['normalization_exact'])==1
assert F(r['U_exact'])==rr['count_U']==7724
s=F(2);tau=F(r['threshold_squared_exact']);cap=F(501,1000)
assert tau==F(1000,499) and s<tau<2*s and 1-s/(2*tau)==cap
assert F(rr['count_threshold_v_exact'])==tau/s
assert F(ri['strict_density_upper'])==F(rr['conservative_expression_upper'])==F('0.024900975')
summary=dict(status='CONDITIONAL',geometry='PASS',radial_primary='PASS',radial_independent='PASS',count='NOT_REPLAYED',dimension=16,cap=str(cap),required_integer_count=7724,strict_conditional_density_upper='0.024900975')
if components_only:
 (D/'COMPONENTS_RESULT.json').write_text(json.dumps(summary,indent=2)+'\n')
 print('CONDITIONAL COMPONENTS PASS; count replay skipped');sys.exit(0)
pair='count/density_kissing16_independent_pair_20261111/'
for f in ['verify.py','replay_chain.py','audit_cyclic.py','replay_gram_integer.py']:
 subprocess.run([sys.executable,str(D/pair/f)],check=True)
c=read(pair+'result.json');chain=read(pair+'replay_result.json');cyc=read(pair+'cyclic_audit_result.json');cover=read(pair+'gram_integer_replay_result.json')
assert c['status'].startswith('PASS') and chain['status']=='PASS' and cyc['status'].startswith('PASS') and cover['status'].startswith('PASS')
assert c['dimension']==16 and F(c['cap_exact'])==cap
assert list(map(F,c['interval_exact']))==[F(-1),cap]
assert F(cover['full_unit_cube_covered_exact'])==1
assert cover['exact_cyclic_bernstein_sha256']==cyc['exact_cyclic_bernstein_sha256']
assert cover['matches_original_partition_transcript'] and cyc['independent_factor_form_matches_exact_cyclic_polynomial']
from check_linkage import count_linkage
objective=F(c['objective_exact']);integer_count=count_linkage(objective,r['U_exact'])
for f,h in c['input_sha256'].items():
 assert h==hashlib.sha256((D/'count/density_kissing16_reweight_frozen_20261111'/f).read_bytes()).hexdigest()
summary.pop('strict_conditional_density_upper',None)
summary.update(status='PASS',count='PASS primary and independent exact replay',integer_count=integer_count,count_objective_exact=str(objective),strict_density_upper='0.024900975',density_interval=rr['density_expression_interval'])
(D/'RESULT.json').write_text(json.dumps(summary,indent=2)+'\n')
print('FULL STRONG16 PASS',json.dumps(summary))
