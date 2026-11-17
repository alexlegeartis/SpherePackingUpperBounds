from pathlib import Path
import json,itertools,numpy as np
from flint import fmpz_mat,fmpq_mat
D=Path(__file__).resolve().parent;z=np.load(D/'lattice80.npz',allow_pickle=False)
G=z['e8_gram'];E=z['e8_twice_euclidean_basis'];F=z['L1_root_basis'];F2=z['L2_root_basis'];T1=z['T1'];T2=z['T2'];C=z['binary_code'];Cp=z['binary_dual'];LC=z['LC_root_basis'];Lp=z['Lprime_root_basis'];Ld=z['Ldoubleprime_twice_root_basis'];M=z['twice_root_coordinate_basis'];GG=z['gram'];Z=z['diagonal_lifts'];Aq=z['old_quotient_basis'];Bq=z['new_quotient_basis']
assert all(np.issubdtype(z[k].dtype,np.integer)for k in z.files)
def det(A):return int(fmpz_mat(A.tolist()).det())
def coords(A,B):
 Q=fmpq_mat(np.array(A).tolist())*fmpq_mat(np.array(B).tolist()).inv();assert all(x.denominator==1 for row in Q.tolist()for x in row);return np.array([[int(x)for x in row]for row in Q.tolist()])
def span(T):return np.array([np.array(a)@T%2 for a in itertools.product(range(2),repeat=len(T))])
assert np.array_equal(E@E.T,4*G)and det(G)==1
assert len(set(map(tuple,span(np.r_[T1,T2]))))==256
for T,Fi,name in [(T1,F,'L1'),(T2,F2,'L2')]:
 coords(T,Fi);coords(2*np.eye(8,dtype=int),Fi);assert abs(det(Fi))==16
 iso=json.loads((D/f'{name}_isometry.json').read_text());S=np.array(iso['basis_change']);Gstd=np.array(iso['canonical_E8_gram']);Est=np.array(iso['canonical_E8_twice_euclidean_basis'])
 assert np.array_equal(Est@Est.T,4*Gstd)and abs(det(S))==1 and np.array_equal(S@Fi@G@Fi.T@S.T,2*Gstd)
 assert not np.any(np.einsum('ij,jk,ik->i',span(T),G,span(T))%4)
assert np.array_equal(C,np.array([[1,0,0,0,0,1,1,1,0,0],[0,1,0,0,0,0,1,1,1,0],[0,0,1,0,0,0,0,1,1,1],[0,0,0,1,0,1,0,0,1,1],[0,0,0,0,1,1,1,0,0,1]]))
for Q in [C,Cp]:assert sorted(span(Q).sum(1))==[0]+[4]*15+[6]*15+[10]
assert not np.any(C@Cp.T%2)
rows=[2*r for r in np.eye(80,dtype=int)]+[np.concatenate([a*t for a in c])for T,Q in [(T1,C),(T2,Cp)]for t in T for c in Q]
coords(rows,LC);assert abs(det(LC))==2**40
coords(Lp,LC);assert abs(det(Lp))==2**44
s=Lp.reshape(80,10,8).sum(1);coords(s,2*F)
# Exact kernel index: induced LC -> E8/L1 has rank4.
s0=LC.reshape(80,10,8).sum(1);assert not np.any(s0%2)
q=s0//2@G@T1.T%2; seen={tuple([0]*4)}
for r in q:seen|={tuple((np.array(x)^r).tolist())for x in list(seen)}
assert len(seen)==16
coords(2*LC,Ld);coords(np.array([np.tile(v,10)for v in F]),Ld);assert abs(det(Ld))==2**116
coords(2*Lp,M);coords(np.array([np.tile(v,10)for v in Z]),M);assert abs(det(M))==2**120
G1=F@G@F.T; assert np.array_equal(Z,Bq@F)
assert len(set(map(tuple,span(np.r_[Aq,Bq]))))==256
old=coords(2*np.eye(8,dtype=int),F)%2
assert set(map(tuple,span(old)))==set(map(tuple,span(Aq)))
for Q in [Aq,Bq]:assert not np.any(np.einsum('ij,jk,ik->i',span(Q),G1,span(Q))%8)
block=np.kron(np.eye(10,dtype=int),G)
assert np.array_equal(M@block@M.T,8*GG)and not np.any(np.diag(GG)%2)and det(GG)==1
assert int(z['minimum8_coefficients']@GG@z['minimum8_coefficients'])==8
for a in itertools.product(range(2),repeat=4):
 v=np.array(a)@Z;coords([2*np.tile(v,10)],2*Lp) # double of half-diagonal lift is in Lprime
 assert tuple(coords([v],F)[0]%2)in set(map(tuple,span(Bq)))
print('PASS exact80D Gram, all intermediate lattices, two E8 isometries, singular complement,16 diagonal shifts, norm8 witness. Minimum≥8 invokes Bachoc–Nebe Theorem4.7.')
