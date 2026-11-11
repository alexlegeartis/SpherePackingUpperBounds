"""Portable independent degree47 upper proof and degree31 floor comparisons."""
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
for script in ['audit_angular_cyclic.py','audit_radial.py','audit_saturation_comparison.py']:
    subprocess.run([sys.executable,str(D/script)],cwd=D,check=True)
results={name:json.loads((D/(name+'_result.json')).read_text())
         for name in ['angular','cyclic','radial','saturation_comparison']}
assert all(r['pass'] for r in results.values())
r=results['radial']
out={'pass':True,'dimension':r['dimension'],'strict_upper':r['strict_density_upper'],
     'density_interval':r['density_interval'],'correction_degree':47,
     'source_t_nonnegative':True,'exact_physical_and_Fourier_span_membership':True,
     'outside_degree31_span_exact':r['outside_degree31_span_exact'],
     'cyclic_rectangles':results['cyclic']['certified_rectangles'],
     'pair_intervals':r['coupling_intervals'],'pair_depth':r['coupling_depth'],
     'all_embedded_lower_proofs_replayed':True,
     'degree31_floors':results['saturation_comparison']['degree31_floors'],
     'free_scalar_amplitude_extension':'SATURATION_EXTENSION.md',
     'input_sha256':pins,'independent_computer_assisted_audit':True,
     'formal_proof_assistant_verification':False}
(D/'audit_result.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='input_sha256'},indent=2))
