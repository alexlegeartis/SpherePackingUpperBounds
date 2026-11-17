from pathlib import Path
import sys,subprocess
n=int(sys.argv[1]);P=Path(__file__).resolve().parent/f'dim{n}'
for name in ['verify_fourier_tail.py','verify_contact_cap.py']+(['verify_angular_structure.py','verify_angular_compact.py']if n==51 else[])+['certify_annuli.py']:
 cmd=[sys.executable,'-B',str(P/name)]
 if name in ['verify_angular_structure.py','verify_angular_compact.py']:cmd.append(str(P))
 subprocess.run(cmd,check=True)
print('FULL REPLAY PASS',n)
