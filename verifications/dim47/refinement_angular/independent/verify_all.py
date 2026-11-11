"""Replay this independent audit from frozen inputs, without an optimizer."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

if not __debug__:
    raise RuntimeError('Do not use Python -O')
D=Path(__file__).resolve().parent
pins={'angular47.json': 'abf4a19d849fd8831b458ab0b56e4a3a09c1bf496067a08435f43a7f627ff17f', 'new_complete_certificate.json': '411ddeccc5425a85e555706e93867b180518318c6c238a47239d87802c85e6cd', 'old_transfer47.json': '116f4457e00b2751295370ad62163f1facadf0a70ed87b73b9d9a186b7b4135d', 'source48_radial.json': '0c87ec795b6d7df3da8e606aa660ed866615a41f4d1b1036f4786c2db4bf977c'}
for name,digest in pins.items():
    assert hashlib.sha256((D/'inputs'/name).read_bytes()).hexdigest()==digest
for script in ['audit_angular_cyclic.py','audit_radial_substitution.py']:
    subprocess.run([sys.executable,str(D/script)],cwd=D,check=True)
results={name:json.loads((D/(name+'_result.json')).read_text())
         for name in ['angular','cyclic','radial_substitution']}
assert all(x['pass'] for x in results.values())
out={'pass':True,'dimension':47,'strict_upper':'0.0000006117584752',
     'density_interval':results['radial_substitution']['density_interval'],
     'cyclic_rectangles':results['cyclic']['certified_rectangles'],
     'pair_intervals':results['radial_substitution']['coupling_intervals'],
     'input_sha256':pins,'independent_computer_assisted_audit':True,
     'formal_proof_assistant_verification':False}
(D/'audit_result.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
