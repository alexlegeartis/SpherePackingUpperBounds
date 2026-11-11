"""Portable radial-only independent replay, with explicit angular dependency."""
from pathlib import Path
import hashlib,json,subprocess,sys
if not __debug__:raise RuntimeError('Do not use Python -O')
D=Path(__file__).resolve().parent
pins=json.loads((D/'INPUT_SHA256.json').read_text())
for name,wanted in pins.items():assert hashlib.sha256((D/name).read_bytes()).hexdigest()==wanted,name
assert pins['inputs/certificate.json']=='25ead1a73f05106192623c3e7cdefa8885cbd607e903523d0ec8716319b4c7b9'
assert pins['inputs/angular_exact.json']=='ebeb0fd1b0532445cf274f01f6f01f1996f608fc9794eb09e3cb43f3932f6354'
subprocess.run([sys.executable,str(D/'audit_radial.py')],cwd=D,check=True)
r=json.loads((D/'radial_result.json').read_text());assert r['pass']
assert not r['full_angular_cyclic_geometry_replayed_here']
out={'pass':True,'scope':'Independent radial audit and PSD/theorem linkage; not a standalone complete angular/cyclic certificate.',
     'dimension':15,'strict_density_expression_upper':r['strict_density_upper'],
     'density_interval':r['density_interval'],'exact_export_and_M_guard':True,
     'independent_Fourier_identity':True,'both_complete_sign_rays':True,
     'pair_intervals':r['pair_intervals'],'pair_depth':r['pair_depth'],
     'pair_domain':r['pair_interval'],'origin_PSD_exact':True,
     'normalization_exact':True,'rational_odd_dimensional_rounding':True,
     'input_sha256':pins,'angular_proof_dependency':'ANGULAR_DEPENDENCY.json and AUDIT.md',
     'full_angular_geometry_replayed_here':False,'formal_proof_assistant_verification':False}
(D/'audit_result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
