"""Fail-closed exact reconstruction and complete minimum replay (C++17 required)."""
from pathlib import Path
import hashlib, subprocess, sys, tempfile
if not __debug__:
    raise RuntimeError('Assertions are required; do not use python -O')
r = Path(__file__).resolve().parent
for line in (r/'SHA256SUMS').read_text().splitlines():
    digest, name = line.split(maxsplit=1)
    if hashlib.sha256((r/name).read_bytes()).hexdigest() != digest:
        raise RuntimeError('Integrity failure: '+name)
subprocess.run([sys.executable, str(r/'generate_and_check.py')], check=True, cwd=r)
with tempfile.TemporaryDirectory() as temp:
    binary = str(Path(temp)/'verify_minimum')
    subprocess.run(['g++', '-std=c++17', '-O2', '-fno-fast-math', '-ffp-contract=off', '-frounding-math', str(r/'interval_enum.cpp'), '-o', binary], check=True)
    result = subprocess.run([binary, str(r/'ldl.txt'), str(r/'Q_reduced.txt')], check=True, text=True, capture_output=True)
    print(result.stdout, end='')
    if result.stdout.strip() != 'COMPLETE 1 nodes 118119521 leaves 1':
        raise RuntimeError('Unexpected or incomplete enumeration: '+result.stdout+result.stderr)
print('PASS: integrity, exact construction and complete minimum proof')
