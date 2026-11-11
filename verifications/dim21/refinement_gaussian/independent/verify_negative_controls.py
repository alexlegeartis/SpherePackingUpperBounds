"""Ensure optimized Python and altered certificate inputs are rejected."""
from pathlib import Path
import json, shutil, subprocess, sys, tempfile
if not __debug__:
    raise RuntimeError('Assertions are required')
root = Path(__file__).resolve().parent
results = {}
for name in ['verify_all.py', 'verify_radial.py', 'verify_angular_cyclic.py', 'verify_comparison.py']:
    p = subprocess.run([sys.executable, '-O', str(root/name)], capture_output=True, text=True)
    assert p.returncode != 0 and 'RuntimeError' in p.stderr, name
    results[f'{name}_rejects_optimized_python'] = True
with tempfile.TemporaryDirectory(prefix='negative_control_', dir=root) as tmp:
    tmp = Path(tmp)
    for name in ['verify_all.py', 'source_hashes.json', 'certificate.json', 'angular.json']:
        shutil.copyfile(root/name, tmp/name)
    with (tmp/'certificate.json').open('a') as f:
        f.write('\n')
    p = subprocess.run([sys.executable, str(tmp/'verify_all.py')], capture_output=True, text=True)
    assert p.returncode != 0 and 'AssertionError: certificate.json' in p.stderr
    results['altered_certificate_rejected_before_proof_execution'] = True
out = {'pass': True, 'dimension': 21, 'checks': results}
(root/'negative_controls.json').write_text(json.dumps(out,indent=2)+'\n')
print(out)
