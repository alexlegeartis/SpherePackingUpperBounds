"""Replay the standalone independent 22D and 23D packing-bound audit."""
import hashlib,json,subprocess,sys
from fractions import Fraction as F
from pathlib import Path
import flint
if not __debug__: raise RuntimeError('Assertions must be enabled')
D=Path(__file__).resolve().parent
assert flint.__version__=='0.9.0'
hashes=json.loads((D/'source_hashes.json').read_text())
source_matches={}
for name,record in hashes.items():
    actual=hashlib.sha256((D/name).read_bytes()).hexdigest()
    assert actual==record['sha256'],name
    source=D.parent.parent/record['source']
    if source.exists():
        assert hashlib.sha256(source.read_bytes()).hexdigest()==actual,str(source)
        source_matches[name]=True
for n in (22,23):
    for verifier in ('verify_angular.py','verify_radial.py'):
        subprocess.run([sys.executable,str(D/verifier),str(n)],cwd=D,check=True)
result={'pass':True,'python_flint':flint.__version__,'source_byte_comparison':source_matches,'bounds':{}}
for n in (22,23):
    radial=json.loads((D/f'radial{n}.json').read_text())
    angular=json.loads((D/f'angular{n}_result.json').read_text())
    comparator=json.loads((D/f'comparators{n}.json').read_text())
    assert radial['pass'] and angular['pass'] and comparator['dimension']==n
    upper=F(radial['strict_upper']); prior=F(comparator['upper']); assert 0<upper<prior
    result['bounds'][n]={'strict_upper':radial['strict_upper'],'pinned_comparator':comparator['upper'],
                         'relative_decrease_lower_exact':str(1-upper/prior),
                         'relative_decrease_lower_percent_approx':float(100*(1-upper/prior))}
(D/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
