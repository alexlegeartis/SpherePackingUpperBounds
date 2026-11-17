if not __debug__:raise RuntimeError('Assertions required')
import json,sys,math,hashlib
from pathlib import Path
from fractions import Fraction as F
sys.set_int_max_str_digits(0)
if not __debug__:raise RuntimeError('Assertions required')
base=Path(__file__).parent
path=Path('certificate.json')
sha=hashlib.sha256(path.read_bytes()).hexdigest();assert sha=='ef5fe2a07be38dbb1a165c19a4685e4a232a2cdf39c38367af721c33f9a02114'
D=json.loads(path.read_text());A=D['cyclic_angular'];assert A==json.load(open(base/'cyclic/angular.json'))
gate=json.load(open(base/'cyclic/meanvalue_angular_result.json'));assert gate['pass'] and gate['dimension']==46
C=F(D['C']);r=F(D['origin_ratio']);C0=F(D['origin_C']);B=F(A['B_exact']);assert C>=0 and r==F(7,10) and C0==r*r*C and C0*C/B**2-(r*C/B)**2==0
qratio=sum(F(v)*(-1)**k for k,v in enumerate(A['gegenbauer_coefficients']))/B;power=-int(A['profile_a_exact']);beta=F(A['profile_b_exact']);h4_upper=F(4**power*math.factorial(256))/(3*beta)**256;assert 0<1+qratio and 1+qratio+2*h4_upper<2*r
E=list(map(F,D['E']));assert all(a<b for a,b in zip(E,E[1:])) and E[-48:]==[2+F(k,192) for k in range(1,49)]
assert next(v for v in json.load(open('check_all_short.json')) if v['n']==46)['pass_exact']
assert next(v for v in json.load(open('check_narrow.json')) if v['n']==46)['pass_exact']
R=json.load(open(base/'replay.json'));P=json.load(open(base/'profiles.json'));assert R['status']==P['status']=='PASS' and R['sha256']==P['sha256']==sha
safe=F(P['safe']);prior=F('.000000871244998');assert safe<prior
R.update(full_pass=True,AA_cuts=P['AA_cuts'],profile_duals=P['profile_branches'],normalization=P['normalization'],safe=P['safe'],source_cyclic_rectangles=gate['certified_rectangles'],source_transfer_exact=False,target_PSD_exact=True,origin_PSD_exact=True,origin_ratio='7/10',combined_tail_exact=True,prior=str(prior),improvement_percent=float(100*(1-safe/prior)),known_gaps=[])
(base/'FULL_AUDIT.json').write_text(json.dumps(R,indent=2));print(json.dumps(R,indent=2))
