"""Aggregate independent proof replay; fresh count replay happens in a temporary copy."""
if not __debug__:raise RuntimeError('Assertions required')
import json,sys,subprocess,shutil,tempfile,hashlib,os
from pathlib import Path
from fractions import Fraction as F
sys.set_int_max_str_digits(0)
D=Path(__file__).resolve().parent;R=D.parent
GEOM=R/'dim16_log_hinge_geometry_independent_20261111';RAD=R/'dim16_log_count_radial_independent_20261111'
COUNT='density_kissing16_reweight_frozen_20261111';PAIR='density_kissing16_independent_pair_20261111';CYCLIC='density_kissing16_independent_cyclic_20261111'
def sh(p):return hashlib.sha256(p.read_bytes()).hexdigest()
expected={
 COUNT+'/candidate.json':'b3225ad4cce819395ce3bd47613f15697af4646fb722ded1916020b5e52ba573',
 COUNT+'/factors.json':'6860eca2c82fdb1b1ebb5d3a4108a216d46c900d6fac7bdeb0d880272e027954',
 CYCLIC+'/polynomial.json':'24085cd076fe5a8d9a4a5bda3b51a199e04b7ea9e82d7e922ce91a0d138bd3e7',
 CYCLIC+'/proof_events.json':'d1d5333e9768677e53f924a0897358567e4164b7a09a6361046072572a669ca5',
 'dim16_log_count_radial_independent_20261111/candidate_exact.json':'3c6d8b01a559b52c4418cccca0e93dde7d34118e4cda4a8f398346f91a22be95'}
for path,hash in expected.items():assert sh(R/path)==hash,path
assert sh(GEOM/'fit.json')==sh(RAD/'inputs/fit.json')=='438c5f0a7d55eca26e6b4c5adbba0c3ec542f0bb18cda15314596456ff2cd4aa'
assert sh(GEOM/'angular.json')==sh(RAD/'inputs/angular_dim16.json')=='e6f6d85ca4fb160536f3ebc74d94f9ab7ea659793cad94258c8a089447a3f40b'

def run(script,log):
 env=os.environ.copy();env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
 with open(D/log,'w')as out:subprocess.run([sys.executable,'-u',str(script)],stdout=out,stderr=subprocess.STDOUT,check=True,env=env)
 print('PASS',script.name,log,flush=True)

with tempfile.TemporaryDirectory(prefix='density16_count_')as td:
 tmp=Path(td)
 for name in [COUNT,PAIR,CYCLIC]:shutil.copytree(R/name,tmp/name,ignore=shutil.ignore_patterns('__pycache__'))
 run(tmp/PAIR/'verify.py','count_pair.log')
 run(tmp/PAIR/'audit_cyclic.py','count_reconstruction.log')
 run(tmp/CYCLIC/'replay_exact.py','count_integer_replay.log')
 for src,dst in [(tmp/PAIR/'result.json','count_pair_replay.json'),(tmp/PAIR/'cyclic_audit_result.json','count_polynomial_reconstruction.json'),(tmp/CYCLIC/'exact_replay_result.json','count_integer_replay.json')]:shutil.copyfile(src,D/dst)

run(GEOM/'verify.py','geometry.log');run(GEOM/'negative_controls.py','geometry_negative.log')
run(RAD/'negative_controls.py','radial_and_negative.log') # Imports and fully runs radial verifier first.
read=lambda p:json.loads(p.read_text())
a=read(D/'count_pair_replay.json');c=read(D/'count_integer_replay.json');g=read(GEOM/'VERIFIED.json');r=read(RAD/'VERIFIED.json');rad=read(RAD/'candidate_exact.json');rnd=read(RAD/'ROUNDING.json')
assert a['dimension']==16 and F(a['cap_exact'])==F(5003,10000)
U=F(a['objective_exact']);integer=U.numerator//U.denominator;assert integer==7802
assert c['verified']and c['events']==2155 and c['covered_unit_cube_volume_exact']=='1'
assert read(D/'count_polynomial_reconstruction.json')['exact_permutation_symmetry']
assert g['status']==r['status']=='PASS'
assert F(rad['U_exact'])==integer and F(rad['threshold_squared_exact'])==F(10000,4997)
s=F(2);tau=F(rad['threshold_squared_exact']);assert s<tau<2*s and 1-s/(2*tau)==F(a['cap_exact'])
assert F(rnd['density_upper_exact'])<F('0.024920988')<F('0.02492121')
neg1=read(GEOM/'NEGATIVE_CONTROLS.json');neg2=read(RAD/'NEGATIVE_CONTROLS.json');assert neg1['status']==neg2['status']=='PASS'
assert all(x['rejected']for x in list(neg1['controls'].values())+list(neg2['controls'].values()))
result={'status':'PASS','dimension':16,'strict_density_upper':'0.024920988','underlying_expression':rnd['density_interval_arb_for_display'],'count_integer':integer,'count_real_exact':str(U),'spherical_cap':'5003/10000','minimum_separation_squared':'2','shell_threshold_squared':str(tau),'geometry_rectangles':g['exact_binary_cover']['accepted'],'geometry_exclusions':g['exact_binary_cover']['discarded'],'count_proof_events':2155,'independent_pair_intervals':r['pair']['leaves'],'radial_rays':'exact rational Sturm,31 variations at both endpoints on each ray','negative_controls':len(neg1['controls'])+len(neg2['controls']),'published_displayed_comparator':'0.02492121','comparator_source':'https://arxiv.org/html/2512.04829v2#S5','proof_input_sha256':expected,'math_scope':'Direct periodic cyclic summation plus pointwise positive-width spherical count; all finite and infinite domains covered; no optimizer status used'}
(D/'COMBINED_RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
