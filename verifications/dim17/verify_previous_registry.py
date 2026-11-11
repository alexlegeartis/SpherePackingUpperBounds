from pathlib import Path
import subprocess,sys,json,hashlib
if not __debug__: raise RuntimeError('Assertions must be enabled')
p=Path(__file__).resolve().parent
assert hashlib.sha256((p/'certificate.json').read_bytes()).hexdigest()=='32a9c22f599c91d80d35a585cbe8ca2e40d53a189d23d6f7ba8c400835f307e4'
r=next(r for r in json.loads((p.parent.parent/'RESULTS.json').read_text()) if r['dimension']==17)
assert r['strict_density_upper_bound']=='0.018335552' and r['package']=='verifications/dim17'
subprocess.run([sys.executable,str(p/'producer_verify.py'),str(p/'certificate.json')],check=True)
subprocess.run([sys.executable,str(p/'independent_verify.py')],check=True)
print('PASS both full-domain cyclic verification paths in 17D')
