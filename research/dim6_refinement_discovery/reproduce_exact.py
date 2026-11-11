"""Rebuild the frozen rational mixture from the saved numerical checkpoint."""
if not __debug__:
    raise RuntimeError('Verification must not run under python -O.')
import argparse
import hashlib
from pathlib import Path
import subprocess
import sys

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('output', type=Path)
args = ap.parse_args()
root = Path(__file__).resolve().parent
out = args.output.resolve()
if out.exists():
    ap.error('Output directory already exists')
out.mkdir(parents=True)
prepared = out / 'continue5_round3_exact.json'
final = out / 'mix_continue_99hundredths.json'
subprocess.run([sys.executable, str(root / 'prepare_candidate.py'),
                str(root / 'continue5_round3_snapshot.json'), str(prepared)], check=True)
subprocess.run([sys.executable, str(root / 'mix_candidates.py'),
                str(root / 'candidate_certificate_v2_6.json'), str(prepared),
                '99/100', str(final)], check=True)
digest = hashlib.sha256(final.read_bytes()).hexdigest()
expected = '41c60fae103355323abdfab8740c64b1089732bc0e0903fd27da845ed4259810'
if digest != expected:
    raise RuntimeError(f'Frozen-source mismatch: {digest}')
print('FROZEN_SOURCE_REPRODUCED', digest)
print('Run export_and_check.py with the proof package for complete verification.')
