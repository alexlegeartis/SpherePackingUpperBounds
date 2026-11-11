"""Replay registered13D primary and independent angular/radial proofs."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,subprocess,sys
if not __debug__:raise RuntimeError('Assertions required')
p=Path(__file__).resolve().parent
r=next(x for x in json.loads((p.parents[1]/'RESULTS.json').read_text())if x['dimension']==13)
assert Fraction(r['strict_density_upper_bound'])==Fraction('0.059905026')and r['package']=='verifications/dim13'
for f,h in json.loads((p/'source_hashes.json').read_text()).items():assert hashlib.sha256((p/f).read_bytes()).hexdigest()==h,f
for f in ['verify_primary.py','independent/verify.py']:subprocess.run([sys.executable,str(p/f)],cwd=p,check=True)
print('PASS registered13D complete primary and independent chains')
