"""Full refined 6D primary and independent proof chains."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,subprocess,sys
if not __debug__:raise RuntimeError('Assertions must remain enabled')
p=Path(__file__).resolve().parent;r=p.parents[2]
row=next(x for x in json.loads((r/'RESULTS.json').read_text()) if x['dimension']==6)
assert Fraction(row['strict_density_upper_bound'])==Fraction('0.409128533') and row['package']=='verifications/dim6/refinement_patch5'
for f,h in json.loads((p/'source_hashes.json').read_text()).items():
 assert hashlib.sha256((p/f).read_bytes()).hexdigest()==h,f
pre=json.loads((p/'prerequisites.json').read_text())
assert Fraction(pre['guards']['reciprocal_gaussian_guards']['eta'])==Fraction(pre['guards']['reciprocal_gaussian_guards']['zeta'])==0
for f in ['verify_primary.py','independent/verify_all.py']:
 subprocess.run([sys.executable,str(p/f)],cwd=p,check=True)
print('PASS improved6D complete primary and independent proof chains, no Gaussian repair')
