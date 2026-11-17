"""Extract the immutable published source package for this numerical experiment."""
from pathlib import Path
import hashlib,tarfile
D=Path(__file__).parent
p=D.parents[1]/'verifications/general52_rank3mid_20261117/proof.tar.xz'
assert hashlib.sha256(p.read_bytes()).hexdigest()=='f8253aa1ad962bd2efcd9c795ba69f646fa829d76896d8e816d110233ff8ab59'
(D/'base').mkdir(exist_ok=True)
with tarfile.open(p) as t:t.extractall(D/'base',filter='data')
print('Verified source archive extracted; this does not certify new candidates.')
