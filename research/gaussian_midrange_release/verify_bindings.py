"""Check frozen evidence and changed integration, not a fresh proof replay."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,subprocess,sys
import flint
if not __debug__:raise RuntimeError('Assertions must remain enabled')
R=Path(__file__).resolve().parents[2];D=Path(__file__).resolve().parent
assert flint.__version__=='0.9.0'
a=json.loads((D/'AUDIT_EQUIVALENCE.json').read_text());rows=json.loads((R/'RESULTS.json').read_text())
for f,h in a['files_sha256'].items():assert hashlib.sha256((R/f).read_bytes()).hexdigest()==h,f
for result in a['results']:
 n=result['dimension'];P=R/f'verifications/dim{n}';row=next(x for x in rows if x['dimension']==n)
 primary=json.loads((P/'full_replay_result.json').read_text());radial=json.loads((P/'independent/audit.json').read_text());angular=json.loads((P/'independent/angular_cyclic_audit.json').read_text())
 assert primary['pass']and radial['pass']and angular['pass']
 assert Fraction(row['strict_density_upper_bound'])==Fraction(primary['strict_upper'])==Fraction(radial['strict_upper_decimal'])
 assert json.loads((R/f'research/dim{n}_gaussian_discovery/inputs/radial_exact.json').read_text())==json.loads((P/'certificate.json').read_text())
 assert json.loads((R/f'research/dim{n}_gaussian_discovery/inputs/angular_exact.json').read_text())==json.loads((P/'continuum_three_point/cyclic_angular40_exact.json').read_text())
for f in [R/'comparison.py',D/'test_parallel_runner.py',D/'test_registration.py']:
 subprocess.run([sys.executable,str(f)],check=True)
print('PASS frozen audit/evidence bindings and changed integration checks; no fresh aggregate mathematical replay is claimed.')
