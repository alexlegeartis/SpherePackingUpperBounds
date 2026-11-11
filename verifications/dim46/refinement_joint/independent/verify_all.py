"""Replay this independent audit from frozen inputs, without an optimizer."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

if not __debug__:
    raise RuntimeError('Do not use Python -O')
D=Path(__file__).resolve().parent
pins={'angular46.json': 'd1577cd3b71b0d80ec39cc1589c66bd73c2c99867d420434334721eb6d55070e', 'certificate.json': '9482f4f13c061f53a38ea3fd8a6668c52d2c91c8a56ff501272942ff2deaeecd', 'discovery_candidate.json': '15a2ee9bf08a0b8764428f716ba7cd3c166c9b59ebecf88ab62e32616f9b3a3d', 'source48_radial.json': '0c87ec795b6d7df3da8e606aa660ed866615a41f4d1b1036f4786c2db4bf977c'}
for name,digest in pins.items():
    assert hashlib.sha256((D/'inputs'/name).read_bytes()).hexdigest()==digest
for script in ['audit_angular_cyclic.py','audit_radial.py']:
    subprocess.run([sys.executable,str(D/script)],cwd=D,check=True)
results={name:json.loads((D/(name+'_result.json')).read_text())
         for name in ['angular','cyclic','radial']}
assert all(x['pass'] for x in results.values())
out={'pass':True,'dimension':46,'strict_upper':'0.0000008868492132',
     'density_interval':results['radial']['density_interval'],
     'cyclic_rectangles':results['cyclic']['certified_rectangles'],
     'pair_intervals':results['radial']['coupling_intervals'],
     'input_sha256':pins,'independent_computer_assisted_audit':True,
     'formal_proof_assistant_verification':False}
(D/'audit_result.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
