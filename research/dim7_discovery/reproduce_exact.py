"""Reproduce the frozen 7D exact source from its saved checkpoint."""
if not __debug__: raise RuntimeError('Verification must not run under python -O.')
import argparse,hashlib,subprocess,sys
from pathlib import Path
root=Path(__file__).resolve().parent
ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('output',type=Path);args=ap.parse_args()
out=args.output.resolve()
if out.exists():ap.error('Output directory already exists')
out.mkdir(parents=True)
prepared=out/'patch_dim7_exact.json'
subprocess.run([sys.executable,str(root/'prepare_candidate.py'),str(root/'candidate_patch_dim7.json'),str(prepared)],check=True)
final=prepared
digest=hashlib.sha256(final.read_bytes()).hexdigest()
if digest!='1cf76290693c7461ad4da8e5389db283ea5bce02c152ea7ab0115bbf88948bf1':raise RuntimeError('Frozen source mismatch: '+digest)
print('FROZEN_SOURCE_REPRODUCED',digest)
