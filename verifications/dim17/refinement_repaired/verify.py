from pathlib import Path
import subprocess,sys,json,hashlib
if not __debug__: raise RuntimeError('Assertions must be enabled')
p=Path(__file__).resolve().parent
assert hashlib.sha256((p/'certificate.json').read_bytes()).hexdigest()=='1bbc96f096ad55b58f2b63276618aa1350798589d09e0e6cf26687cef1900405'
r=next(r for r in json.loads((p.parent.parent.parent/'RESULTS.json').read_text()) if r['dimension']==17)
assert r['strict_density_upper_bound']=='0.018204096' and r['package']=='verifications/dim17/refinement_repaired'
subprocess.run([sys.executable,str(p/'verify_repaired.py'),str(p/'certificate.json')],check=True)
subprocess.run([sys.executable,str(p/'independent_verify.py')],check=True)
print('PASS mixed-Gaussian source, cyclic geometry and exact Schur saturation: 17D <0.018204096')
