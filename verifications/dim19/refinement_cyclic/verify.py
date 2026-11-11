from pathlib import Path
import subprocess,sys,json,hashlib
if not __debug__: raise RuntimeError('Assertions must be enabled')
p=Path(__file__).resolve().parent
assert hashlib.sha256((p/'certificate.json').read_bytes()).hexdigest()=='9a125b3b65c1845e2857e19004b2731b76f0b963096433a968a79f4575900db3'
r=next(r for r in json.loads((p.parent.parent.parent/'RESULTS.json').read_text()) if r['dimension']==19)
assert r['strict_density_upper_bound']=='0.009557491' and r['package']=='verifications/dim19/refinement_cyclic'
subprocess.run([sys.executable,str(p/'producer_verify.py'),str(p/'certificate.json')],check=True)
subprocess.run([sys.executable,str(p/'independent_verify.py')],check=True)
print('PASS both full-domain cyclic verification paths in 19D')
