"""Complete43D independently implemented exact/interval proof audit."""
from pathlib import Path
import hashlib,json,subprocess,sys
if not __debug__:raise RuntimeError('Assertions required')
p=Path(__file__).resolve().parent
for f,h in json.loads((p/'source_hashes.json').read_text()).items():
 assert hashlib.sha256((p/f).read_bytes()).hexdigest()==h,f
for f in ['verify_angular_cyclic.py','verify_radial.py']:
 subprocess.run([sys.executable,str(p/f)],check=True,cwd=p)
print('PASS complete43D exact angular/cyclic/radial/origin/rounding audit')
