from pathlib import Path
import numpy as np,json,hashlib,itertools,random
from flint import fmpz_mat as Z
D=Path(__file__).parent;v=list(map(int,(D/'input.txt').read_text().split()));pos=0
size=v[pos];pos+=1;stored=np.array(v[pos:pos+size*8],dtype=np.int64).reshape(size,8);pos+=size*8;cs=v[pos:pos+size];pos+=size;J=np.array(v[pos:pos+64],dtype=np.int64).reshape(8,8)
assert size==9121 and len(set(map(tuple,stored)))==size
actual=set()
def rec(parity,arr,total):
 if len(arr)==8:
  if sum(arr)%4==0:actual.add(tuple(arr))
  return
 for t in range(-4,5):
  if t%2==parity and total+t*t<=24:rec(parity,arr+[t],total+t*t)
rec(0,[],0);rec(1,[],0)
assert actual==set(map(tuple,stored))
counts={k:sum(sum(t*t for t in x)==4*k for x in actual)for k in [0,2,4,6]};assert counts=={0:1,2:240,4:2160,6:6720}
# Explicit standard-coordinate E8 basis and index, independently of producer Gram claims.
source=np.load(D.parent/'bachoc_nebe80_reconstruction_20261117/lattice80.npz');E=source['e8_twice_euclidean_basis'];assert abs(int(Z(E.tolist()).det()))==256
for row in E:assert len({int(x)%2 for x in row})==1 and int(sum(row))%4==0
assert np.array_equal(E@J,4*np.eye(8,dtype=np.int64))
# Exhaustive 512 coordinate candidates against a full shell scan on60 small targets.
rng=random.Random(811701);checks=[]
for a in [[0]*8,[8]*8,[-8]*8]+[[rng.randrange(-8,9)for _ in range(8)]for _ in range(57)]:
 best=10**9
 for off in [0,8]:
  lower=[(x-off)//16 for x in a]
  for word in range(256):
   zz=[lower[i]+((word>>i)&1)for i in range(8)]
   if sum(zz)%2==0:best=min(best,sum((a[i]-16*zz[i]-off)**2 for i in range(8)))
 delta=np.array(a)-8*stored;direct=int(np.min(np.sum(delta*delta,axis=1)));assert best==direct and best<=256;checks.append(best)
# Exact local radius implication: (1+sqrt(5/2))²<8 since sqrt(10)<9/2.
assert 10<81/4
out=dict(status='PASS independent complete shell, E8 coordinate basis and nearest-point controls',counts=counts,controls=len(checks),input_sha=hashlib.sha256((D/'input.txt').read_bytes()).hexdigest())
(D/'local_result.json').write_text(json.dumps(out,indent=2));print(out)
