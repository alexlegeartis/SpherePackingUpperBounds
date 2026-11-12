if not __debug__:raise RuntimeError('Assertions required')
from pathlib import Path
from fractions import Fraction as F
import subprocess,sys,json
D=Path(__file__).resolve().parent
dep=json.loads((D/'external_dependencies.json').read_text());assert dep['dimension']==9 and dep['average_ordered_neighbor_bound']==373
assert F(dep['used_radius_ratio'])==F('1.003012')<F(dep['published_radius_ratio'])==F('1.003013')
assert dep['source']=='https://arxiv.org/html/1206.2608#S6.T2' and dep['original_sdp_replayed_here'] is False
for name in ['fit.json','angular.json']:
 for folder in ['radial','independent_radial','independent_geometry']:
  assert (D/'geometry'/name).read_bytes()==(D/folder/name).read_bytes(),('LINKAGE',folder,name)
for folder in ['radial','independent_radial']:
 assert (D/'certificate.json').read_bytes()==(D/folder/'certificate.json').read_bytes()
assert F(json.loads((D/'certificate.json').read_text())['threshold_v_exact'])==F(dep['used_radius_ratio'])**2
for script in ['geometry/verify_geometry.py','radial/verify.py','independent_radial/verify_radial_independent.py','independent_geometry/verify_geometry_independent.py','linkage/verify_linkage.py']:
 subprocess.run([sys.executable,str(D/script)],check=True)
for name in ['geometry/RESULT_GEOMETRY.json','radial/RESULT.json','independent_radial/RADIAL_INDEPENDENT.json','independent_geometry/VERIFIED.json']:
 assert json.loads((D/name).read_text())['status']=='PASS'
print('PASS: density9 <0.190926892; published average-count373 external')
