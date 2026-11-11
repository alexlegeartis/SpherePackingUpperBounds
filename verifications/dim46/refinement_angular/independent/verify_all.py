"""Replay this independent audit from frozen inputs, without an optimizer."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

if not __debug__:
    raise RuntimeError('Do not use Python -O')
D=Path(__file__).resolve().parent
pins={
    'angular46.json':'d1577cd3b71b0d80ec39cc1589c66bd73c2c99867d420434334721eb6d55070e',
    'radial46.json':'d26fb1e751f1a475d1955531ef04b8e120a1827299707bb37fd1f4ae1aeb1aaa',
    'new_complete_certificate.json':'75861548c8d1376075aa3cf324ad6760e4b4b56caad8938a9eca1d2ca6c6d6d6',
}
for name,digest in pins.items():
    assert hashlib.sha256((D/'inputs'/name).read_bytes()).hexdigest()==digest
for script in ['audit_angular_cyclic.py','audit_radial_substitution.py']:
    subprocess.run([sys.executable,str(D/script)],cwd=D,check=True)
results={name:json.loads((D/(name+'_result.json')).read_text())
         for name in ['angular','cyclic','radial_substitution']}
assert all(x['pass'] for x in results.values())
out={'pass':True,'dimension':46,'strict_upper':'0.000000897267583',
     'density_interval':results['radial_substitution']['density_interval'],
     'cyclic_rectangles':results['cyclic']['certified_rectangles'],
     'pair_intervals':results['radial_substitution']['coupling_intervals'],
     'input_sha256':pins,'independent_computer_assisted_audit':True,
     'formal_proof_assistant_verification':False}
(D/'audit_result.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
