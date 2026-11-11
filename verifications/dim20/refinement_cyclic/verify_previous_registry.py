from pathlib import Path
import subprocess,sys,json,hashlib
if not __debug__: raise RuntimeError('Assertions must be enabled')
p=Path(__file__).resolve().parent
assert hashlib.sha256((p/'certificate.json').read_bytes()).hexdigest()=='69647964d2d0ecc93b47cada16809ea879dc9561811463140d76c5d8407ead3e'
r=next(r for r in json.loads((p.parent.parent.parent/'RESULTS.json').read_text()) if r['dimension']==20)
assert r['strict_density_upper_bound']=='0.006977702' and r['package']=='verifications/dim20/refinement_cyclic'
subprocess.run([sys.executable,str(p/'producer_verify.py'),str(p/'certificate.json')],check=True)
subprocess.run([sys.executable,str(p/'independent_verify.py')],check=True)
print('PASS both full-domain cyclic verification paths in 20D')
