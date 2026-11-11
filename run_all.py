"""Replay every registered proof package. Never run with python -O."""
import subprocess,sys
from pathlib import Path
if not __debug__: raise RuntimeError('Assertions are required')
r=Path(__file__).resolve().parent
subprocess.run([sys.executable,str(r/'comparison.py')],check=True)
jobs=[('dim6','verify.py'),('dim45','verify.py'),('dim46','verify.py'),('dim47','verify.py'),('dim17','refinement_repaired/verify.py'),('dim19','refinement_cyclic/verify.py'),('dim20','refinement_cyclic/verify.py'),('dim17','verify.py'),('dim18','refinement_cyclic/verify.py'),('dim18','verify.py'),('dim5','refinement_p4/verify.py'),('dim4','refinement_p6/verify.py'),('dim5','refinement_pd2/verify.py'),('dim4','refinement_p4/verify.py'),('dim5','verify.py'),('dim4','verify.py'),('dim21','verify.py'),('dim21','independent_verify.py'),('dim19','verify.py'),('dim19','independent_verify.py'),('dim20','verify.py'),('dim20','independent_verify.py')]+[('dim48','continuum_three_point/independent_cyclic40_'+x+'.py') for x in ['angular','domain','radial']]
for directory,script in jobs:
 print(directory,script,flush=True)
 subprocess.run([sys.executable,script],cwd=r/'verifications'/directory,check=True)
from fractions import Fraction
import json
a=json.loads((r/'verifications/dim48/continuum_three_point/independent_cyclic40_radial.json').read_text())
assert Fraction(a['density_upper']) < Fraction('4.13385739314e-7')
subprocess.run([sys.executable,str(r/'verifications/origin_slack/verify_algebra.py')],check=True)
subprocess.run([sys.executable,str(r/'verifications/cyclic_origin_slack/verify_algebra.py')],check=True)
print('PASS: all registered proof replays and strict decimal rounding')
