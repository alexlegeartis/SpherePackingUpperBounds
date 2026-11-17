from pathlib import Path
import sys,subprocess,argparse
D=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--quick',action='store_true');args=ap.parse_args()
def run(path,*argv):subprocess.run([sys.executable,'-B',str(D/path),*argv],check=True)
if not args.quick:
 run('replay_trace.py','main');run('replay_trace.py','cap');run('input/root_repair/certify_annular.py');run('input/root_repair/check_tail.py')
run('check_exact.py')
print('PASS lattice48 exact algebra'+(' with archived interval receipts'if args.quick else ' and complete geometry/radial interval replay'))
