"""Fail-closed controls run against temporary copies; frozen inputs untouched."""
if not __debug__:raise RuntimeError('Assertions required')
from pathlib import Path
from fractions import Fraction as F
import tempfile,shutil,json,subprocess,sys
D=Path(__file__).resolve().parent
results={}
with tempfile.TemporaryDirectory(prefix='log_hinge_controls_')as td:
 T=Path(td)
 def run(name,mutate):
  for fn in ['verify.py','fit.json','angular.json','partition.json','RESULT_GEOMETRY.json']:shutil.copyfile(D/fn,T/fn)
  mutate(T)
  r=subprocess.run([sys.executable,str(T/'verify.py')],capture_output=True,text=True,timeout=60)
  assert r.returncode!=0 and ('AssertionError'in r.stderr or 'RuntimeError'in r.stderr),name
  results[name]={'rejected':True}
 def missing(t):
  p=json.loads((t/'partition.json').read_text());p['accepted'].pop();(t/'partition.json').write_text(json.dumps(p))
 def false_exclusion(t):
  p=json.loads((t/'partition.json').read_text());j=next(i for i,z in enumerate(p['accepted'])if 2*F(z['box'][0])*F(z['box'][2])<=1);z=p['accepted'].pop(j);p['discarded'].append(z['box']);(t/'partition.json').write_text(json.dumps(p))
 def negative_profile(t):
  p=json.loads((t/'fit.json').read_text());p['parameters'][7]=-1;(t/'fit.json').write_text(json.dumps(p))
 run('missing_geometry_rectangle',missing);run('invalid_sorted_domain_exclusion',false_exclusion);run('negative_profile_coefficient',negative_profile)
 r=subprocess.run([sys.executable,'-O',str(D/'verify.py')],capture_output=True,text=True)
 assert r.returncode!=0 and 'Assertions required'in r.stderr
 results['python_optimization_disabled']={'rejected':True}
(D/'NEGATIVE_CONTROLS.json').write_text(json.dumps({'status':'PASS','controls':results},indent=2)+'\n');print(results)
