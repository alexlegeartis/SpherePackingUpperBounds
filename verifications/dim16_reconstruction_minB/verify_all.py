"""Full preserved base plus reconstructed 16D transfer proof."""
if not __debug__: raise RuntimeError('Assertions required')
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,subprocess,sys
D=Path(__file__).resolve().parent
row=next(r for r in json.loads((D.parents[1]/'RESULTS.json').read_text()) if r['dimension']==16)
assert F(row['strict_density_upper_bound'])==F('0.024920883') and row['package']=='verifications/dim16_reconstruction_minB'
for f,h in json.loads((D/'input_hashes.json').read_text()).items():
 assert hashlib.sha256((D/f).read_bytes()).hexdigest()==h,('HASH',f)
for f in ['fit.json','angular.json']:
 assert (D/f).read_bytes()==(D/'independent/geometry'/f).read_bytes(),('LINKAGE',f)
for f in ['../dim16/verify.py','verify_geometry.py','verify_scaled_bound.py','independent/geometry/verify.py','independent/verify_transfer.py']:
 subprocess.run([sys.executable,str(D/f)],check=True)
a=json.loads((D/'TRANSFER_RESULT.json').read_text());b=json.loads((D/'independent/TRANSFER_AUDIT.json').read_text())
assert F(a['strict_upper'])==F(b['bound'])==F('0.024920883')
assert F(a['C_exact'])==F(b['C_exact'])
print('PASS: full primary and independent reconstructed 16D transfer, including published base replays')
