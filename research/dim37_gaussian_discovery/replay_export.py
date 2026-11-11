"""Deterministic export replay, not an optimizer run or a substitute for verify.py."""
from pathlib import Path
import subprocess,sys,json
D=Path(__file__).resolve().parent
subprocess.run([sys.executable,str(D/'export_angular.py'),'37'],check=True)
assert json.loads((D/'work/cyclic_angular40_exact.json').read_text())==json.loads((D/'inputs/angular_exact.json').read_text())
subprocess.run([sys.executable,str(D/'export_radial.py'),str(D/'inputs/radial_scaled37_degree31_r4p625.json')],check=True)
got=json.loads((D/'work/radial_scaled37_degree31_r4p625_certificate.json').read_text());want=json.loads((D/'inputs/radial_exact.json').read_text());assert got==want
print('Frozen 37D angular and radial exports exactly reproduced; no optimizer rerun claimed.')
