"""Replay inherited continuum certificates, then the exact Schur tightening."""
from pathlib import Path
import subprocess, sys
assert __debug__, 'Assertions are required'
P=Path(__file__).resolve().parent
import json, hashlib
for name,digest in json.loads((P/'source_inputs.json').read_text()).items():
 assert hashlib.sha256((P.parent.parent/name).read_bytes()).hexdigest()==digest, name
jobs=[('dim4/refinement_p6','verify.py'),('dim5/refinement_p4','verify.py')]+[(f'dim{n}',s) for n in [19,20,21] for s in ['verify.py','independent_verify.py']]
for directory,script in jobs:
 print('Source replay:',directory,script,flush=True)
 subprocess.run([sys.executable,script],cwd=P.parent/directory,check=True)
subprocess.run([sys.executable,str(P/'verify_algebra.py')],check=True)
print('PASS all source domains, Schur factorization and strict outward rounding')
