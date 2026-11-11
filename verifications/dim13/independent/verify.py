"""Replay both independent mathematical checks on pinned inputs."""
from pathlib import Path
import hashlib,json,subprocess,sys
if not __debug__:raise RuntimeError('Do not use -O or -OO')
D=Path(__file__).resolve().parent
pins=json.loads((D/'pinned_inputs.json').read_text())
for name,digest in pins.items():assert hashlib.sha256((D/'inputs'/name).read_bytes()).hexdigest()==digest,name
for script in ['audit_angular.py','audit_radial.py']:
 subprocess.run([sys.executable,str(D/script)],cwd=D,check=True)
a=json.loads((D/'independent_result.json').read_text());r=json.loads((D/'radial_result.json').read_text())
assert a['pass'] and r['pass']
for name,digest in pins.items():assert hashlib.sha256((D/'inputs'/name).read_bytes()).hexdigest()==digest,name
print('COMPLETE_INDEPENDENT_AUDIT_PASS',r['dimension'],r['strict_density_upper'])
