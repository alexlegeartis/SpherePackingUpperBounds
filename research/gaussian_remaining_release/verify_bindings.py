"""Check saved proof bindings and changed integration behavior, not a fresh proof replay."""
from pathlib import Path
import hashlib,json,subprocess,sys
if not __debug__:raise RuntimeError('Assertions required')
p=Path(__file__).resolve().parent;r=p.parents[1]
for name,key in [('NEW_PROOF_HASHES.json',None),('HISTORICAL_FILES.json','sha256')]:
 data=json.loads((p/name).read_text());data=data[key]if key else data
 for f,h in data.items():assert hashlib.sha256((r/f).read_bytes()).hexdigest()==h,f
results=json.loads((p/'PRIMARY_EXPORT_RESULTS.json').read_text())
assert len(results)==6 and all(x['primary_pass'] and x['export_pass'] for x in results)
results=json.loads((p/'INDEPENDENT_RESULTS.json').read_text())
assert len(results)==6 and all(x['independent_full_chain_pass'] for x in results)
for f in [r/'comparison.py',p/'test_registration.py',p/'test_parallel_runner.py']:
 subprocess.run([sys.executable,str(f)],check=True)
print('PASS recorded proof hashes and changed integration checks; not a fresh full proof replay')
