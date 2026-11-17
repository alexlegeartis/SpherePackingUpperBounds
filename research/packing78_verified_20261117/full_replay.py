from pathlib import Path
import sys,subprocess
P=Path(__file__).resolve().parent
for f in ['mw80_normal162_20261117/parent/verify_parent.py','packing78_parent_independent_20261117/check_parent.py','packing78_independent_20261117/check.py','packing78_independent_20261117/check_external.py']:
 subprocess.run([sys.executable,'-B',str(P/f)],check=True)
print('PASS complete78D independent parent, geometry, density and external example')
