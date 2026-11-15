from pathlib import Path
import subprocess,sys,json,hashlib
if not __debug__: raise RuntimeError('Assertions must be enabled')
p=Path(__file__).resolve().parent
assert hashlib.sha256((p/'certificate.json').read_bytes()).hexdigest()=='1bbc96f096ad55b58f2b63276618aa1350798589d09e0e6cf26687cef1900405'
r=next(r for r in json.loads((p.parent.parent.parent/'RESULTS.json').read_text()) if r['dimension']==17)
from fractions import Fraction
if r['package']=='verifications/dim17/refinement_repaired':
 assert Fraction(r['strict_density_upper_bound'])==Fraction('0.018204096')
else:
 assert r['package']=='verifications/dim17/refinement_larger20261115'
 assert 'verifications/dim17/refinement_repaired' in r.get('earlier_repository_packages',[])
 assert Fraction(r['strict_density_upper_bound'])<Fraction('0.018204096')
subprocess.run([sys.executable,str(p/'verify_repaired.py'),str(p/'certificate.json')],check=True)
subprocess.run([sys.executable,str(p/'independent_verify.py')],check=True)
print('PASS mixed-Gaussian source, cyclic geometry and exact Schur saturation: 17D <0.018204096')
