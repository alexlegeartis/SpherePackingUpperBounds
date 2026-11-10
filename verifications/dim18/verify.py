from pathlib import Path
import subprocess,sys
if not __debug__: raise RuntimeError('Assertions must be enabled')
p=Path(__file__).resolve().parent
import hashlib,json
assert hashlib.sha256((p/'certificate.json').read_bytes()).hexdigest()=='286d466bfd83e3ec8543f7fa7312655432701289a51f84c76a369afbdf0ee137'
r=next(r for r in json.loads((p.parent.parent/'RESULTS.json').read_text()) if r['dimension']==18)
assert r['strict_density_upper_bound']=='0.013408931' and r['package']=='verifications/dim18'
subprocess.run([sys.executable,str(p/'run_certificate.py'),str(p/'certificate.json')],check=True)
subprocess.run([sys.executable,str(p/'independent_verify.py')],check=True)
print('PASS: both 18D verification paths')
