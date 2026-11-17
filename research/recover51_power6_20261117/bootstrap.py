"""Extract the immutable published source package for this numerical experiment."""
from pathlib import Path
import hashlib,tarfile
D=Path(__file__).parent
p=D.parents[1]/'verifications/general51_rank3p8_20261117/proof.tar.xz'
assert hashlib.sha256(p.read_bytes()).hexdigest()=='7d1b48627948429301d916f3ed04beb255d33c7cd325cb47a7ca81acbb7128e6'
(D/'base').mkdir(exist_ok=True)
with tarfile.open(p) as t:t.extractall(D/'base',filter='data')
print('Verified source archive extracted; this does not certify new candidates.')
