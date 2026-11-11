from pathlib import Path
import subprocess,sys,json,hashlib
if not __debug__: raise RuntimeError('Assertions must be enabled')
p=Path(__file__).resolve().parent
assert hashlib.sha256((p/'certificate.json').read_bytes()).hexdigest()=='84fe1f8779e84c4bb7f4643ad7afc9afcc9761dd9f793fd0f33e934319055b82'
r=next(r for r in json.loads((p.parent.parent.parent/'RESULTS.json').read_text()) if r['dimension']==18)
from fractions import Fraction
assert Fraction(r['strict_density_upper_bound']) <= Fraction('0.013284728')
assert r['package'] in ['verifications/dim18/refinement_cyclic','verifications/cyclic_origin_slack']
subprocess.run([sys.executable,str(p/'producer_verify.py'),str(p/'certificate.json')],check=True)
subprocess.run([sys.executable,str(p/'independent_verify.py')],check=True)
print('PASS both full-domain cyclic verification paths in 18D')
