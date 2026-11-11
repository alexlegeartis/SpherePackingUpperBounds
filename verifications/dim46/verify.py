"""Replay inherited angular/cyclic proof, new radial proof and independent audit."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,subprocess,sys
if not __debug__:raise RuntimeError('Assertions must remain enabled')
p=Path(__file__).resolve().parent
r=p.parents[1]
rows=json.loads((r/'RESULTS.json').read_text());row=next(x for x in rows if x['dimension']==46)
assert Fraction(row['strict_density_upper_bound'])<=Fraction('0.000000914882847')
assert row['package'] in ['verifications/dim46','verifications/dim46/refinement_angular']
for name,h in json.loads((p/'source_hashes.json').read_text()).items():
 assert hashlib.sha256((p/name).read_bytes()).hexdigest()==h,name
subprocess.run([sys.executable,str(p/'verify_primary.py')],check=True)
subprocess.run([sys.executable,str(p/'independent/verify_all.py')],check=True)
print('PASS registered46D complete primary and independent proof chains')
