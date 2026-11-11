"""Reproduce the frozen 5D exact source from its saved checkpoint."""
if not __debug__: raise RuntimeError('Verification must not run under python -O.')
import argparse,hashlib,subprocess,sys
from pathlib import Path
root=Path(__file__).resolve().parent
ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('output',type=Path);args=ap.parse_args()
out=args.output.resolve()
if out.exists():ap.error('Output directory already exists')
out.mkdir(parents=True)
prepared=out/'patch_dim5_exact.json'
subprocess.run([sys.executable,str(root/'prepare_candidate.py'),str(root/'candidate_patch_dim5.json'),str(prepared)],check=True)
final=out/'patch_dim5_guarded_1e8.json'
subprocess.run([sys.executable,str(root/'add_kernel_guard.py'),str(prepared),str(root/'kernel_guard5.json'),'1/100000000',str(final)],check=True)
digest=hashlib.sha256(final.read_bytes()).hexdigest()
if digest!='79881d16bc8d837d13222c0d3c5ca6d47264e031e09cd46738275198d4da85b1':raise RuntimeError('Frozen source mismatch: '+digest)
print('FROZEN_SOURCE_REPRODUCED',digest)
