"""Complete registered 6D source and independent audit replay."""
from pathlib import Path
from fractions import Fraction
import json,hashlib,subprocess,sys
if not __debug__:raise RuntimeError('Assertions must remain enabled')
p=Path(__file__).resolve().parent
row=next(x for x in json.loads((p.parents[1]/'RESULTS.json').read_text()) if x['dimension']==6)
assert Fraction(row['strict_density_upper_bound'])==Fraction('0.410189987') and row['package']=='verifications/dim6'
for f,h in json.loads((p/'source_hashes.json').read_text()).items():
 assert hashlib.sha256((p/f).read_bytes()).hexdigest()==h,f
for f in ['verify_primary.py', 'audit_independent.py', 'independent/verify_algebra_guards.py', 'independent/verify_pair.py']:
 subprocess.run([sys.executable,str(p/f)]+([str(p)] if f=='audit_independent.py' else []),cwd=p,check=True)
print('PASS registered 6D full producer and independent audit chains')
