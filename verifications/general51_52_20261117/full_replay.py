from pathlib import Path
import sys,subprocess,argparse
P=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--quick',action='store_true');a=ap.parse_args()
if not a.quick:
 for n in [51,52]:subprocess.run([sys.executable,'-B',str(P/'run_dimension.py'),str(n)],check=True)
for name in ['check_contact.py','check_exact.py']:subprocess.run([sys.executable,'-B',str(P/name)],check=True)
print('PASS exact51/52 checks'+(' using saved interval receipts'if a.quick else ' and fresh full interval replay'))
