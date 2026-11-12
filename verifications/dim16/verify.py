"""Replay registered16D primary and independent angular/radial proofs."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,subprocess,sys
if not __debug__:raise RuntimeError('Assertions required')
p=Path(__file__).resolve().parent
r=next(x for x in json.loads((p.parents[1]/'RESULTS.json').read_text())if x['dimension']==16)
assert Fraction('0.024920988') in [Fraction(r['strict_density_upper_bound'])]+[Fraction(x) for x in r.get('earlier_repository_upper_densities',[])], 'Preserved base bound absent from registry history'
for f,h in json.loads((p/'source_hashes.json').read_text()).items():assert hashlib.sha256((p/f).read_bytes()).hexdigest()==h,f
for f in ['verify_all.py','independent/dim16_log_count_combined_independent_20261111/run_all.py']:subprocess.run([sys.executable,str(p/f)],cwd=p,check=True)
print('PASS preserved16D base complete primary and independent chains')
