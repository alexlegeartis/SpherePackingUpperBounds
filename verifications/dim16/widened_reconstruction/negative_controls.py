"""Corrupt disposable copies only; verify rejection."""
if not __debug__:raise RuntimeError('Assertions required')
from pathlib import Path
from tempfile import TemporaryDirectory
import shutil,subprocess,sys,json
from check_linkage import count_linkage
D=Path(__file__).resolve().parent;results={}
for name,f in [('count_hash','count/density_kissing16_reweight_frozen_20261111/candidate.json'),('radial_coefficient','radial/candidate_exact.json')]:
 with TemporaryDirectory() as td:
  target=Path(td)/'package';shutil.copytree(D,target)
  p=target/f;j=json.loads(p.read_text())
  if name=='radial_coefficient':j['physical'][0]='0'
  else:j['corrupted_for_negative_control']=True
  p.write_text(json.dumps(j))
  r=subprocess.run([sys.executable,str(target/'verify_all.py')],capture_output=True,text=True)
  assert r.returncode!=0 and 'HASH' in r.stderr,(name,r.stdout,r.stderr)
  results[name]='REJECTED by pinned mathematical input hash'
c=json.loads((D/'count/density_kissing16_independent_pair_20261111/result.json').read_text())
assert count_linkage(c['objective_exact'],'7724')==7724
try:count_linkage(c['objective_exact'],'7723')
except AssertionError:results['understated_count7723']='REJECTED by exact floor/radial linkage'
else:raise AssertionError('Invalid count accepted')
(D/'NEGATIVE_CONTROLS.json').write_text(json.dumps(results,indent=2)+'\n');print(results)
