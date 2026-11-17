from pathlib import Path
import sys,subprocess,argparse
D=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--quick',action='store_true');ap.add_argument('--geometry',required=True);args=ap.parse_args()
def run(path,*argv):subprocess.run([sys.executable,'-B',str(path),*argv],check=True)
run(Path(args.geometry)/'run_all.py',*(['--quick']if args.quick else []))
if not args.quick:
 for n in [46,47]:
  run(D/f'dim{n}/certify_annular.py');run(D/f'dim{n}/check_tail.py')
run(D/'check_exact.py')
print('PASS compact46/47 transfer'+(' using saved interval receipts'if args.quick else ' and full inherited geometry plus scalar replays'))
