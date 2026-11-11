"""Replay the registered 46D primary and independent proof chains."""
from pathlib import Path
from fractions import Fraction
import json,hashlib,subprocess,sys
if not __debug__:raise RuntimeError('Assertions must remain enabled')
p=Path(__file__).resolve().parent
row=next(x for x in json.loads((p.parents[2]/'RESULTS.json').read_text()) if x['dimension']==46)
assert Fraction(row['strict_density_upper_bound'])<=Fraction('0.000000897267583') and row['package'] in ['verifications/dim46/refinement_angular','verifications/dim46/refinement_joint','verifications/dim46/refinement_degree47','verifications/dim46/refinement_extended_pair']
for f,h in json.loads((p/'source_hashes.json').read_text()).items():
 assert hashlib.sha256((p/f).read_bytes()).hexdigest()==h,f
for f in ['verify_primary.py','independent/verify_all.py']:
 subprocess.run([sys.executable,str(p/f)],cwd=p,check=True)
print('PASS registered 46D complete primary and independent chains')
