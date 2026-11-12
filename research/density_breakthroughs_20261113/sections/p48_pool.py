"""Explicit norm-six pool by integral root exchanges, no completeness claim."""
from pathlib import Path
import re,json,sys
import numpy as np
P=Path(__file__).resolve().parent;ROOT=P.parents[2]
EXT=ROOT/'research/density_breakthroughs_20261113/seed_deformations/external_sphere_packing'
ambient=sys.argv[1] if len(sys.argv)>1 else 'P48p'
limit=int(sys.argv[2])if len(sys.argv)>2 else 8192
def matrix_after(text,key):
 s=text.index(key);a=text.index('[[',s);b=text.index(']]',a)+2
 return np.array([[int(z)for z in re.findall(r'-?\d+',line)]for line in text[a:b].splitlines()],dtype=np.int64)
paths=sorted((EXT/'packings_out').glob(f'**/{ambient}-*.cert.txt'))
assert paths
G=matrix_after(paths[0].read_text(),'ambient Gram G')
assert G.shape==(48,48)and np.array_equal(G,G.T)
vectors=[];seen=set()
def add(v):
 nz=np.flatnonzero(v)
 if not len(nz):return
 if v[nz[0]]<0:v=-v
 key=tuple(map(int,v))
 if key not in seen:
  assert int(v@G@v)==6
  seen.add(key);vectors.append(v.copy())
for v in np.eye(48,dtype=np.int64):
 if v@G@v==6:add(v)
for path in paths:
 text=path.read_text()
 assert np.array_equal(matrix_after(text,'ambient Gram G'),G)
 T=matrix_after(text,'coefficient rows w.r.t.')
 for v in T:
  if v@G@v==6:add(v)
rng=np.random.default_rng(130036)
for it in range(10000):
 n=len(vectors)
 if n>=limit:break
 V=np.array(vectors)
 # Many triples x,y,x-y lie in the norm-six shell when <x,y>=3.
 a=rng.integers(n,size=20000);b=rng.integers(n,size=20000)
 inn=np.einsum('bi,ij,bj->b',V[a],G,V[b]);ok=np.where(np.abs(inn)==3)[0]
 for j in ok:
  add(V[a[j]]-int(np.sign(inn[j]))*V[b[j]])
  if len(vectors)>=limit:break
 if it%10==0:print(ambient,it,len(vectors),flush=True)
 if len(vectors)==n and it>20:break
V=np.array(vectors,dtype=np.int64)
assert np.all(np.einsum('bi,ij,bj->b',V,G,V)==6)
np.savez_compressed(P/f'{ambient}_pool.npz',G=G,V=V)
(P/f'{ambient}_pool.json').write_text(json.dumps({'ambient':ambient,'source_certificates':[str(x.relative_to(ROOT))for x in paths],'vectors_up_to_sign':len(V),'exact_norm':6,'complete':False,'basis':'ambient Gram from cert; all vectors integer row coefficients'},indent=2))
print('DONE',ambient,len(V),flush=True)
