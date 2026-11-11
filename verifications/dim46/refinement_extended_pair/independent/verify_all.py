"""Portable independent 46D extended-pair proof. No optimizer or network."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

if not __debug__:
    raise RuntimeError('Do not use Python -O')
D=Path(__file__).resolve().parent
pins=json.loads((D/'INPUT_SHA256.json').read_text())
for name,wanted in pins.items():
    assert hashlib.sha256((D/name).read_bytes()).hexdigest()==wanted,name
assert pins['inputs/certificate.json']=='1eb9563cadf892101a24b87e5e9bfb2686b1c216297ebe29521ae32f506c6a01'
assert pins['inputs/angular46.json']=='d1577cd3b71b0d80ec39cc1589c66bd73c2c99867d420434334721eb6d55070e'
for script in ['audit_angular_cyclic.py','audit_radial.py']:
    subprocess.run([sys.executable,str(D/script)],cwd=D,check=True)
results={name:json.loads((D/(name+'_result.json')).read_text())
         for name in ['angular','cyclic','radial','tail_escape']}
assert all(r['pass'] for r in results.values())
r=results['radial'];e=results['tail_escape']
assert r['coupled_upper_X_exact']==r['physical_tail_X_exact']=='40'
assert r['M_guard_exact']=='10001/10000'
assert e['old_negative_tail_template_violated_exact']
out={'pass':True,'dimension':46,'strict_upper':r['strict_density_upper'],
     'density_interval':r['density_interval'],'correction_degree':47,
     'source_t_nonnegative':True,'exact_physical_and_Fourier_span_membership':True,
     'positive_H0_normalization_checked':True,'M_guard_exact':r['M_guard_exact'],
     'complete_Fourier_ray':True,'physical_tail_X_exact':'40',
     'full_pair_interval':'[6*pi,40]','pair_intervals':r['coupling_intervals'],
     'pair_depth':r['coupling_depth'],'pair_interval_union_checked':True,
     'angular_PSD_and_origin_Gram':True,
     'cyclic_rectangles':results['cyclic']['certified_rectangles'],
     'complete_cyclic_domain_and_tails':True,'exact_rational_density_rounding':True,
     'tail_escape_point_X_exact':e['point_X_exact'],
     'tail_escape_positive_exact':True,'old_negative_tail_template_excluded':True,
     'tail_scope_theorem':'AUDIT.md',
     'input_sha256':pins,'independent_computer_assisted_audit':True,
     'formal_proof_assistant_verification':False}
(D/'audit_result.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='input_sha256'},indent=2))
