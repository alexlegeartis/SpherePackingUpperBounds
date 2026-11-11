"""Standalone radial and JOINT saturation replay; no optimizer needed."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import subprocess
import sys

if not __debug__:
    raise RuntimeError('Do not use Python -O')
D=Path(__file__).resolve().parent
pins=json.loads((D/'INPUT_SHA256.json').read_text())
assert pins['inputs/angular46_lower_certificate.json']=='4f345f28a5738aafb6dca3c5fc13b8e25715440a81eead84bab01c99488f6b45'
assert pins['inputs/angular47_lower_certificate.json']=='6b96cfe7ef6ac626a7d8c566c4feed400d17ca5049d89835df52ca36e94a8a79'
assert pins['inputs/angular46_caponly_certificate.json']=='99b65a54ed50c3bef2634ac37f7297cd1a2903a529e1dccc788b878bcdf93b9a'
assert pins['inputs/angular47_caponly_certificate.json']=='22ea2ccda17c2d121cfecab1d59d6f00a8688993f09d5e281fa8490935f43d0b'
assert pins['verify_angular_certificate.py']=='802385efbd83ff5c5437e2318f7a59f56cfe85ac88c9debc60411c65838ea71c'
for name,wanted in pins.items():
    assert hashlib.sha256((D/name).read_bytes()).hexdigest()==wanted,name

results=[]
for n in [46,47]:
    angular_out=D/f'angular{n}_replayed.json'
    subprocess.run([sys.executable,str(D/'verify_angular_certificate.py'),
                    str(D/f'inputs/angular{n}_lower_certificate.json'),'--out',str(angular_out)],check=True)
    angular=json.loads(angular_out.read_text())
    assert angular['pass'] and angular['all_harmonic_degrees_certified']
    cap_path=D/f'inputs/angular{n}_caponly_certificate.json'
    cap_source=json.loads(cap_path.read_text())
    assert all(row['point'][0]=='cap' for row in cap_source['dual_rows'])
    cap_out=D/f'angular{n}_caponly_replayed.json'
    subprocess.run([sys.executable,str(D/'verify_angular_certificate.py'),str(cap_path),'--out',str(cap_out)],check=True)
    cap=json.loads(cap_out.read_text())
    assert cap['pass'] and cap['all_harmonic_degrees_certified']
    for tag in ['fixed','joint','allrate']:
        subprocess.run([sys.executable,str(D/'verify_dual.py'),str(D/f'dual{n}_{tag}_selected.json')],check=True)
        r=json.loads((D/f'dual{n}_{tag}_verified.json').read_text())
        assert r['pass']
        if tag=='joint':assert F(r['B_exact'])==F(angular['lower_bound_B_exact'])
        if tag=='allrate':
            assert F(r['B_exact'])==F(cap['lower_bound_B_exact'])
            assert r['all_pair_rows_at_contact']
        upper=F('0.0000008868492132' if n==46 else '0.0000006010083731')
        lower=F(r['safe_density_lower_exact'])
        gap=1-lower/upper
        results.append({'dimension':n,'scope':tag,'B_exact':r['B_exact'],
                        'certified_family_density_lower':r['safe_density_lower_decimal'],
                        'comparison_upper':str(upper),'remaining_fractional_improvement_upper_exact':str(gap),
                        'remaining_improvement_ppm_upper':float(gap*10**6),
                        'radial_dual_pass':True,'all_degree_angular_replay_pass':tag!='fixed',
                        'every_positive_scalar_Gaussian_rate':tag=='allrate'})
out={'pass':True,'results':results,'joint_bridge':'B>=L_B and M>=0 imply f0+BM>=f0+L_B*M; radial feasible cone in(f,M) is independent ofB',
     'radial_restrictions':'degree111 source plus Laguerre degrees0..31, t>=0, M>=0, H(0)=1; physical P<=0 beyond7pi; pair coupling on squared radii[6,7]',
     'angular_restrictions':'scalar separable Gaussian kernel; all finite harmonic degrees; alpha=51pi/4 for joint and every positive alpha for allrate',
     'scope_warning':'These are LOWER bounds on the best certificate objective in the specified family, not lower bounds on packing density or unrestricted-method saturation.',
     'comparison_upper_source':'Previously independently verified dim46/47_border_refit packages; this replay verifies the lower certificates, not those separate upper packages.'}
(D/'saturation_result.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
