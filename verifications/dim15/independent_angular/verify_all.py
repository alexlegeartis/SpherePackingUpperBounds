"""Standalone angular/theorem replay; no optimizer or radial inputs required."""
if not __debug__:
    raise RuntimeError('Assertions are required; do not run this verifier with -O or -OO.')

from pathlib import Path
import subprocess,sys,json,hashlib
D=Path(__file__).resolve().parent
subprocess.run([sys.executable,str(D/'independent_check.py')],check=True)
r=json.loads((D/'independent_result.json').read_text());assert r['pass']
out={'pass':True,'angular_input_sha256':r['angular_sha256'],'checker_sha256':hashlib.sha256((D/'independent_check.py').read_bytes()).hexdigest(),'profile':'h(u)=exp[-12(u/s-1)-16(u/s-1)^2], s=15/8','dimension':15,'covers_all_squared_triangle_sides_at_least_s':True,'PSD_angular_and_rank_one_origin_linkage':True,'cap_and_single_sign_transition':True,'full_cyclic_geometry_and_dilation_lemma':True,'radial_certificate_audited_here':False,'density_claim_made_here':False,'details':'independent_result.json; AUDIT.md'}
(D/'AUDIT_RESULT.json').write_text(json.dumps(out,indent=2)+'\n');print(out,flush=True)
