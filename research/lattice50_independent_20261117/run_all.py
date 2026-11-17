from pathlib import Path
import subprocess,sys
P=Path(__file__).resolve().parent
for s in ['replay_annuli.py','producer/verify_tail_transform.py','check_exact.py']:
 subprocess.run([sys.executable,'-B',str(P/s)],check=True)
print('PASS reviewed50D lattice certificate')
