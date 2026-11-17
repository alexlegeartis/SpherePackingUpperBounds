from pathlib import Path
import subprocess,sys,argparse
P=Path(__file__).resolve().parent
if not __debug__:raise RuntimeError('Assertions required')
a=argparse.ArgumentParser();a.add_argument('--quick',action='store_true');args=a.parse_args()
for n in [49,50]:
 D=P/f'dim{n}'
 if not args.quick:
  for name in ['verify_fourier_tail.py','verify_angular_structure.py','verify_angular_compact.py','certify_annuli.py']:
   cmd=[sys.executable,'-B',str(D/name)]
   if name in ['verify_angular_structure.py','verify_angular_compact.py']:cmd.append(str(D))
   subprocess.run(cmd,check=True)
for name in ['check_exact49.py','dim50/check_exact_independent.py','check_caps.py']:
 subprocess.run([sys.executable,'-B',str(P/name)],check=True)
print('PASS49/50 exact algebra and caps'+(' (recorded interval proof dependencies; quick mode)' if args.quick else ' and complete fresh interval proofs'))
