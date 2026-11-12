"""Replay the complete 10D certificate and enforce component linkage."""
if not __debug__:raise RuntimeError('Assertions required')
import subprocess,sys,json,hashlib
from pathlib import Path
D=Path(__file__).parent
for f,h in json.loads((D/'input_hashes.json').read_text()).items():
 assert hashlib.sha256((D/f).read_bytes()).hexdigest()==h,('HASH',f)
for name in ['fit.json','angular.json']:
 assert (D/'geometry'/name).read_bytes()==(D/'radial'/name).read_bytes(),('LINKAGE',name)
 assert (D/'geometry'/name).read_bytes()==(D/'independent'/name).read_bytes(),('INDEPENDENT_LINKAGE',name)
assert (D/'radial/certificate.json').read_bytes()==(D/'independent/certificate.json').read_bytes()
for script in ['geometry/verify_geometry.py','radial/verify.py','independent/verify_geometry_independent.py','independent/verify_radial_independent.py','negative_controls.py']:
 subprocess.run([sys.executable,str(D/script)],check=True)
g=json.loads((D/'geometry/RESULT_GEOMETRY.json').read_text());r=json.loads((D/'radial/RESULT.json').read_text())
assert g['status']==r['status']=='PASS'
assert r['dimension']==10 and r['strict_bound']=='0.143330282'
print('PASS: complete reconstructed 10D certificate, density < 0.143330282')
