"""Replay source48, primary transfer and independent exact transfer."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,subprocess,sys
if not __debug__:raise RuntimeError('Assertions must remain enabled')
p=Path(__file__).resolve().parent
r=p.parents[1]
rows=json.loads((r/'RESULTS.json').read_text());row=next(x for x in rows if x['dimension']==47)
assert Fraction(row['strict_density_upper_bound'])<=Fraction('0.000000616879891')
assert row['package'] in ['verifications/dim47','verifications/dim47/refinement_angular','verifications/dim47/refinement_joint','verifications/dim47/refinement_degree47']
for name,h in json.loads((p/'source_hashes.json').read_text()).items():
 assert hashlib.sha256((p/name).read_bytes()).hexdigest()==h,name
subprocess.run([sys.executable,str(p/'verify_all.py')],check=True)
subprocess.run([sys.executable,str(p/'independent/verify.py')],check=True)
print('PASS registered47D complete primary and independent proof chains')
