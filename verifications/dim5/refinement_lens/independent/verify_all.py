"""Complete independent audit chain for the frozen 5D tail-funded lens witness."""
from pathlib import Path
import subprocess,sys,json,hashlib
if not __debug__:raise RuntimeError('Assertions required')
p=Path(__file__).resolve().parent
for name,digest in json.loads((p/'source_hashes.json').read_text()).items():
 assert hashlib.sha256((p/name).read_bytes()).hexdigest()==digest,name
for script in ['verify_algebra_guards.py','verify_pair.py','verify_budgets.py','certify_cyclic_bernstein.py','verify_far.py']:
 args=[sys.executable,str(p/script)]
 if script=='certify_cyclic_bernstein.py':args.append(str(p/'witness_psd_exact.json'))
 subprocess.run(args,cwd=p,check=True)
print('PASS complete independently audited tail-funded lens density <0.509996147')
