"""Replay registered 23D primary and independent complete proof chains."""
from pathlib import Path
from fractions import Fraction
import json,hashlib,subprocess,sys
if not __debug__:raise RuntimeError('Assertions must remain enabled')
p=Path(__file__).resolve().parent
row=next(x for x in json.loads((p.parents[1]/'RESULTS.json').read_text()) if x['dimension']==23)
assert Fraction(row['strict_density_upper_bound'])==Fraction('0.002642540089') and row['package']=='verifications/dim23'
for f,h in json.loads((p/'source_hashes.json').read_text()).items():
 assert hashlib.sha256((p/f).read_bytes()).hexdigest()==h,f
subprocess.run([sys.executable,str(p/'verify_primary.py')],cwd=p,check=True)
for f in ['verify_angular.py','verify_radial.py']:
 subprocess.run([sys.executable,str(p.parent/'dim22_23_independent'/f),'23'],cwd=p.parent/'dim22_23_independent',check=True)
print('PASS registered 23D complete primary and independent proof chains')
