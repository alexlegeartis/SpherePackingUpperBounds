"""Replay the registered 7D primary and independent full proof chains."""
from pathlib import Path
from fractions import Fraction
import json,hashlib,subprocess,sys
if not __debug__:raise RuntimeError('Assertions must remain enabled')
p=Path(__file__).resolve().parent;r=p.parents[1]
row=next(x for x in json.loads((r/'RESULTS.json').read_text()) if x['dimension']==7)
assert Fraction(row['strict_density_upper_bound'])==Fraction('0.320584772') and row['package']=='verifications/dim7'
for f,h in json.loads((p/'source_hashes.json').read_text()).items():
 assert hashlib.sha256((p/f).read_bytes()).hexdigest()==h,f
for f,args in [('verify_primary.py',[]),('audit_independent.py',[str(p)]),('independent/verify_all.py',[])]:
 subprocess.run([sys.executable,str(p/f)]+args,cwd=p,check=True)
print('PASS registered 7D full primary and independent audit chains')
