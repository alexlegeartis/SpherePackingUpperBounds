"""Registered5D spherical-lens proof, with separate full independent audit."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,subprocess,sys
if not __debug__:raise RuntimeError('Assertions must remain enabled')
p=Path(__file__).resolve().parent;r=p.parents[2]
row=next(x for x in json.loads((r/'RESULTS.json').read_text()) if x['dimension']==5)
assert Fraction(row['strict_density_upper_bound'])==Fraction('0.509953498') and row['package']=='verifications/dim5/refinement_spherical_lens'
for f,h in json.loads((p/'source_hashes.json').read_text()).items():
 assert hashlib.sha256((p/f).read_bytes()).hexdigest()==h,f
assert (p/'source_exact.json').read_bytes()==(p/'independent/source_exact.json').read_bytes()
assert json.loads((p/'step_allowances_spherical_005.json').read_text())==json.loads((p/'independent/budget.json').read_text())
for f in ['verify_primary.py','independent/verify_all.py']:
 subprocess.run([sys.executable,str(p/f)],cwd=p,check=True)
print('PASS registered5D spherical lens full producer and independent chains')
