from pathlib import Path
import json,hashlib,itertools,numpy as np
from flint import fmpq_mat
D=Path(__file__).resolve().parent;P=D.parent/'l80_code_decoder_20261117';L=D.parent/'bachoc_nebe80_reconstruction_20261117';z=np.load(L/'lattice80.npz',allow_pickle=False);G=z['gram'];M=z['twice_root_coordinate_basis'];E=z['e8_twice_euclidean_basis'];F=z['L1_root_basis'];H=np.r_[z['T1'],z['T2']];Z=z['diagonal_lifts'];c=json.loads((P/'case101.json').read_text());r=np.array(c['r'],dtype=np.int64)
assert r.shape==(80,)and set(r.tolist())<={0,1}and np.any(r);assert np.array_equal(np.array(c['gram']),G);assert int(r@G@r)%4==0
v=list(map(int,(P/'case101.txt').read_text().split()));cursor=0
def get(n,shape=None):
 global cursor
 a=np.array(v[cursor:cursor+n],dtype=np.int64);assert len(a)==n;cursor+=n;return a.reshape(shape)if shape else a
n=int(get(1)[0]);assert n==9121;DS=get(8*n,(n,8));DC=get(n);J=get(64,(8,8));lookup=get(256);X=get(80,(10,8));zs=get(128,(16,8));states=get(2048,(256,8));masks=get(25,(5,5));parity=get(32);assert cursor==len(v)
assert np.array_equal(E@J,4*np.eye(8,dtype=int))
lookup_expected={};rootstates=[]
for s in range(256):
 coeff=np.array([(s>>k)&1 for k in range(8)])@H;rootstates.append(coeff);key=sum(int(coeff[j]%2)<<j for j in range(8));assert key not in lookup_expected;lookup_expected[key]=s>>4
assert np.array_equal(lookup,[lookup_expected[i]for i in range(256)])
rootstates=np.array(rootstates);assert np.array_equal(states,rootstates@E)
assert np.array_equal(X,(r@M).reshape(10,8)@E)
zz=np.array([np.array([(s>>k)&1 for k in range(4)])@Z for s in range(16)]);assert np.array_equal(zs,zz@E)
assert np.array_equal(np.array(c['target_euclidean_times8']),X)and np.array_equal(np.array(c['diagonal_z_root']),zz)
# All16 shifts are exact members of the selected complement; kappa=0.
Bq=z['new_quotient_basis'];old=z['old_quotient_basis'];spanB={tuple((np.array(b)@Bq)%2)for b in itertools.product(range(2),repeat=4)}
def ic(A,B):
 q=fmpq_mat(np.array(A).tolist())*fmpq_mat(np.array(B).tolist()).inv();assert all(t.denominator==1 for row in q.tolist()for t in row);return np.array([[int(t)for t in row]for row in q.tolist()])
for root in zz:
 assert tuple(ic([5*root],F)[0]%2)in spanB
 ic([2*np.tile(root,10)],2*z['Lprime_root_basis'])
assert len(set(tuple(root%2)for root in zz))==16
C=z['binary_code'];words=np.array(list(itertools.product(range(2),repeat=10)));dual=words[np.all(words@C.T%2==0,axis=1)];Q=[]
for i in range(5):
 found=dual[np.all(dual[:,:5]==np.eye(5,dtype=int)[i],axis=1)];assert len(found)==1;Q.append(found[0])
Q=np.array(Q);assert np.array_equal(np.array(c['code_C']),C)and np.array_equal(np.array(c['code_dual']),Q)
expected_masks=np.zeros((5,5),dtype=int)
for i in range(5):
 for j in range(5):
  expected_masks[i,j]=sum(int(C[i,5+j])<<b for b in range(4))+sum(int(Q[i,5+j])<<b for b in range(4,8))
assert np.array_equal(masks,expected_masks)
for i in range(32):
 b=np.array([(i>>j)&1 for j in range(5)]);word=b@Q%2;assert int(word.sum())%2==0;assert parity[i]==int(word.sum())//2%2
# Carries of every supplied shell vector are checked independently in root coordinates.
num=DS@J;assert not np.any(num%4)
for coeff,carry in zip(num//4,DC):
 key=sum(int(coeff[j]%2)<<j for j in range(8));assert carry==lookup_expected[key]
# The quadratic target is an integral-distance half-class: distances² are
# rGr/4 + yGy − rGy, integral for any lattice coefficient vector y.
assert int(r@G@r)%4==0
out={'status':'PASS exact case101 input maps, code masks, carries,16glues and target normalization','target_parent_norm':int(r@G@r),'target_squared_norm':int(r@G@r)//4,'local_cost_to_global_norm_denominator':128,'absence_cutoff_numerator':640,'absence_scope':'No absence proof in this checker; cutoff auditor separately checks local/global search completeness.','hashes':{str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in [P/'case101.txt',P/'case101.json',L/'lattice80.npz']}}
(D/'result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
