from pathlib import Path
import numpy as np,itertools,json
from flint import fmpz_mat,fmpq_mat
D=Path(__file__).resolve().parent;z=np.load(D/'lattice80.npz');E=z['e8_twice_euclidean_basis'];F=z['L1_root_basis'];G=z['e8_gram'];T=z['T1']
R=[]
for i in range(8):
 for s in [-4,4]:v=np.zeros(8,dtype=int);v[i]=s;R.append(v)
for inds in itertools.combinations(range(8),4):
 for signs in itertools.product([-2,2],repeat=4):v=np.zeros(8,dtype=int);v[list(inds)]=signs;R.append(v)
for i in range(8):
 for signs in itertools.product([-1,1],repeat=8):
  v=np.array(signs);v[i]*=3
  if sum(v)%4==0:R.append(v)
assert len(R)==2160
Q=fmpq_mat(np.array(R).tolist())*fmpq_mat(E.tolist()).inv();assert all(x.denominator==1 for row in Q.tolist()for x in row)
X=np.array([[int(x)for x in row]for row in Q.tolist()]);sel=np.flatnonzero(np.all(X@G@T.T%2==0,axis=1));assert len(sel)==240
R=np.array(R)[sel];X=X[sel]
positive=R@(10**np.arange(8))>0;pos=R[positive];posset=set(map(tuple,pos));simple=[]
for r in pos:
 if not any(tuple(r-a)in posset for a in pos):simple.append(r)
assert len(simple)==8
Q=fmpq_mat(np.array(simple).tolist())*fmpq_mat(E.tolist()).inv()*fmpq_mat(F.tolist()).inv();assert all(x.denominator==1 for row in Q.tolist()for x in row)
S=np.array([[int(x)for x in row]for row in Q.tolist()]);assert abs(int(fmpz_mat(S.tolist()).det()))==1
J=S@F@G@F.T@S.T//2;assert np.array_equal(np.diag(J),[2]*8)and set(J.flatten())=={-1,0,2}
Estd=np.zeros((8,8),dtype=int);Estd[0]=[1,-1,-1,-1,-1,-1,-1,1];Estd[1,:2]=[2,2]
for j in range(2,8):Estd[j,j-2]=-2;Estd[j,j-1]=2
Gstd=Estd@Estd.T//4
perm=next(p for p in itertools.permutations(range(8))if np.array_equal(J[np.ix_(p,p)],Gstd))
S=S[list(perm)];assert np.array_equal(S@F@G@F.T@S.T,2*Gstd)
(D/'L1_isometry.json').write_text(json.dumps({'basis_change':S.tolist(),'canonical_E8_gram':Gstd.tolist(),'canonical_E8_twice_euclidean_basis':Estd.tolist(),'norm4_roots_in_L1':240,'basis_change_determinant':int(fmpz_mat(S.tolist()).det())},indent=2)+'\n');print('PASS explicit unimodular L1 basis change to twice canonical E8 Gram')
