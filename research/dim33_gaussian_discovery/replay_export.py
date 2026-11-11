"""Deterministic export replay, not an optimizer run or a substitute for verify.py."""
from pathlib import Path
import subprocess,sys,json
D=Path(__file__).resolve().parent
subprocess.run([sys.executable,str(D/'export_angular.py'),'33'],check=True)
assert json.loads((D/'work/cyclic_angular40_exact.json').read_text())==json.loads((D/'inputs/angular_exact.json').read_text())
subprocess.run([sys.executable,str(D/'export_radial.py'),str(D/'inputs/radial_scaled33_degree31_r4p125.json')],check=True)
got=json.loads((D/'work/radial_scaled33_degree31_r4p125_certificate.json').read_text());want=json.loads((D/'inputs/radial_exact.json').read_text());assert got==want
print('Frozen 33D angular and radial exports exactly reproduced; no optimizer rerun claimed.')
