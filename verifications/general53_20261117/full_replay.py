from pathlib import Path
import sys,subprocess,argparse
D=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--quick',action='store_true');args=ap.parse_args()
def run(path,*argv):subprocess.run([sys.executable,'-B',str(D/path),*argv],check=True)
if not args.quick:
 run('dim53/verify_fourier_tail.py');run('dim53/certify_annuli.py','53')
run('check_contact.py');run('check_exact.py')
print('PASS53 exact algebra'+(' with archived interval receipts'if args.quick else ' and full1024-bit interval replay'))
