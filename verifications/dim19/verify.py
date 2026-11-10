from pathlib import Path
import subprocess,sys
root=Path(__file__).resolve().parent
for name in ('verify_global.py','verify_annular.py'):
 subprocess.run([sys.executable,str(root/name),str(root/'certificate.json')],check=True)
print('ALL CERTIFICATE CHECKS PASSED. Density < 0.009618351.')
