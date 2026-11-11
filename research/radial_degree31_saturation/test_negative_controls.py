"""Adversarial validation of row-domain and nonsingularity guards."""
from pathlib import Path
import copy
import json
import subprocess
import sys
import tempfile

if not __debug__:
    raise RuntimeError('Do not use Python -O')
D=Path(__file__).resolve().parent
original=json.loads((D/'dual46_fixed_selected.json').read_text())
cases=[]
def add(name,mutator):
    record=copy.deepcopy(original);mutator(record);cases.append((name,record))
add('unknown_row_kind',lambda r:r['rows'][0].update(kind='unknown'))
add('pair_below_separation',lambda r:next(s for s in r['rows'] if s['kind']=='pair').update(u='5'))
add('negative_multiplier',lambda r:r['rows'][0].update(multiplier_exact='-1'))
def physical_zero(r):
    row=next(s for s in r['rows'] if s['kind']=='physical');row.pop('pi_interpolation',None);row['X']='0'
add('physical_point_outside_tail',physical_zero)
add('negative_fourier_radius',lambda r:next(s for s in r['rows'] if s['kind']=='fourier').update(X='-1'))
add('singular_duplicate_row',lambda r:r['rows'].__setitem__(1,copy.deepcopy(r['rows'][0])))
add('nonpositive_B',lambda r:r.update(B_exact='0'))
results=[]
with tempfile.TemporaryDirectory(dir=D) as tmp:
    for name,r in cases:
        path=Path(tmp)/(name+'.json');path.write_text(json.dumps(r))
        out=subprocess.run([sys.executable,str(D/'verify_dual.py'),str(path)],capture_output=True,text=True)
        assert out.returncode!=0,name
        results.append({'test':name,'rejected':True,'last_error':out.stderr.strip().splitlines()[-1]})
out=subprocess.run([sys.executable,'-O',str(D/'verify_dual.py'),str(D/'dual46_fixed_selected.json')],capture_output=True,text=True)
assert out.returncode!=0
results.append({'test':'optimized_python_disallowed','rejected':True,'last_error':out.stderr.strip().splitlines()[-1]})
(D/'negative_controls.json').write_text(json.dumps({'pass':True,'tests':results},indent=2)+'\n')
print(json.dumps({'pass':True,'rejected_tests':len(results)},indent=2))
