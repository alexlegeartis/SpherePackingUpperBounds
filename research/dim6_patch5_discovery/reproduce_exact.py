"""Reproduce the separately frozen direct follow-on exact source."""
if not __debug__:
    raise RuntimeError('Verification must not run under python -O.')
import argparse
import hashlib
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parent
ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('output', type=Path)
args = ap.parse_args()
out = args.output.resolve()
if out.exists():
    ap.error('Output directory already exists')
out.mkdir(parents=True)
source = out / 'patch5_exact.json'
subprocess.run([sys.executable, str(root / 'prepare_candidate.py'),
                str(root / 'candidate_patch5.json'), str(source)], check=True)
digest = hashlib.sha256(source.read_bytes()).hexdigest()
expected = '19bddaa444bb214613549d2ee1ec6c786c585c4bc1c72f7fcd3977acef4ef521'
if digest != expected:
    raise RuntimeError(f'Frozen-source mismatch: {digest}')
print('FROZEN_SOURCE_REPRODUCED', digest)
