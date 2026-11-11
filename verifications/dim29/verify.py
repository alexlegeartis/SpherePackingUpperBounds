"""Replay the registered 29D primary and independent proof chains."""
from pathlib import Path
from fractions import Fraction
import json,hashlib,subprocess,sys
if not __debug__:raise RuntimeError('Assertions must remain enabled')
p=Path(__file__).resolve().parent
row=next(x for x in json.loads((p.parents[1]/'RESULTS.json').read_text()) if x['dimension']==29)
assert Fraction(row['strict_density_upper_bound'])==Fraction('0.0003420173821') and row['package']=='verifications/dim29'
for f,h in json.loads((p/'source_hashes.json').read_text()).items():
 assert hashlib.sha256((p/f).read_bytes()).hexdigest()==h,f
for f in ['verify_primary.py','independent/verify_all.py']:
 subprocess.run([sys.executable,str(p/f)],cwd=p,check=True)
print('PASS registered 29D complete primary and independent chains')
