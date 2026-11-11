"""Rebuild the frozen exact 5D lens source from its final numerical checkpoint."""
if not __debug__:
    raise RuntimeError('Verification must not run under python -O.')

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
SOURCE_SHA256 = '28585450c46b4fbdc6f9eefec55b7e6d29f18fc527b654d00359ed425cb1f949'
CHECKPOINT_SHA256 = '9b373c621c8a85e99b592cb78102bc09d88c1c6d87a754d53ae2b9a65bd2c5d1'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path, help='New output directory')
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error('Output directory already exists; choose a new path')
    checkpoint = ROOT / 'wedge5_spherical.json'
    if hashlib.sha256(checkpoint.read_bytes()).hexdigest() != CHECKPOINT_SHA256:
        raise RuntimeError('Final numerical checkpoint hash mismatch')
    if hashlib.sha256((ROOT / 'source_exact.json').read_bytes()).hexdigest() != SOURCE_SHA256:
        raise RuntimeError('Bundled frozen source hash mismatch')
    output.mkdir(parents=True)
    exact = output / 'source_exact.json'
    subprocess.run([sys.executable, str(ROOT / 'prepare_candidate.py'),
                    str(checkpoint), str(exact)], check=True)
    digest = hashlib.sha256(exact.read_bytes()).hexdigest()
    if digest != SOURCE_SHA256 or exact.read_bytes() != (ROOT / 'source_exact.json').read_bytes():
        raise RuntimeError('Frozen exact source mismatch: ' + digest)
    record = json.loads(exact.read_text())[0]
    if any(r['off_origin_diagonal_shift'] != '0'
           for r in record['exact_source_preparation']['repairs']):
        raise RuntimeError('An off-origin shift occurred')
    print('FROZEN_SOURCE_REPRODUCED', digest)
    print('BYTE_IDENTICAL; all 13 off-origin shifts are zero. This is not certification.')


if __name__ == '__main__':
    main()
