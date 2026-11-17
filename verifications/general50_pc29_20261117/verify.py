from pathlib import Path
import json,hashlib,subprocess,sys
if not __debug__:raise RuntimeError('Do not use Python -O')
p=Path(__file__).resolve().parent
for name,sha in json.loads((p/'MANIFEST.json').read_text()).items():
 assert hashlib.sha256((p/name).read_bytes()).hexdigest()==sha,name
jobs=[['check_bindings.py'],['geometry/verify_caps.py','50'],['scalar/verify_fourier_tail.py'],['scalar/certify_annuli.py'] if '--full' in sys.argv else None,['scalar/verify_caps.py','50']]
if '--full' in sys.argv:jobs.append(['replay_geometry.py','50'])
for job in jobs:
 if job:subprocess.run([sys.executable,str(p/job[0]),*job[1:]],check=True)
print('PASS certificate checks'+(' and full domain/radial replays'if '--full'in sys.argv else '; use --full to replay interval coverage'))
