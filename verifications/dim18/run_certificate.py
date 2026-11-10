"""Run both independent condition groups; never promote discovery automatically."""
import argparse,subprocess,sys
from pathlib import Path
if not __debug__:raise RuntimeError('Assertions must be enabled')
a=argparse.ArgumentParser();a.add_argument('witness',type=Path);args=a.parse_args()
witness=args.witness.resolve();root=Path(__file__).resolve().parent
for script in ['verify_global.py','verify_annular.py']:
 subprocess.run([sys.executable,str(root/script),str(witness)],check=True)
print('ALL COMPUTATIONAL CONDITIONS PASS. Literature comparison and independent review remain separate.')
