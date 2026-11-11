"""Full replay of diagnostic certificates, not record registration."""
from pathlib import Path
import hashlib,json,subprocess,sys
if not __debug__: raise RuntimeError('Assertions must remain enabled')
p=Path(__file__).resolve().parent
for f,h in json.loads((p/'source_hashes.json').read_text()).items():
    assert hashlib.sha256((p/f).read_bytes()).hexdigest()==h, f
for f in ['run_range_applicability.py','independent16/verify.py','independent22/verify.py','independent22/schur.py']:
    subprocess.run([sys.executable,str(p/f)],check=True)
print('PASS: both primary and independent complete proofs, including origin saturation; NONCOMPETITIVE diagnostics only.')
