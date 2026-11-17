from pathlib import Path
import numpy as np,json,itertools,hashlib
from flint import fmpz_mat,fmpq_mat
D=Path(__file__).resolve().parent
# Twice the standard simple-root coordinates of E8.
E=np.zeros((8,8),dtype=np.int64);E[0]=[1,-1,-1,-1,-1,-1,-1,1];E[1,:2]=[2,2]
for j in range(2,8):E[j,j-2]=-2;E[j,j-1]=2
alignment=json.loads((D/'decoder_alignment.json').read_text()); E=np.array(alignment['E8_basis_times2'],dtype=np.int64)
G=E@E.T//4;assert np.all(E@E.T%4==0)and int(fmpz_mat(G.tolist()).det())==1
V=np.array(list(itertools.product(range(2),repeat=8)),dtype=np.int64)
def span(rows):
 out={tuple([0]*8)}
 for r in rows:out|={tuple((np.array(x)^r).tolist())for x in list(out)}
 return out
def independent(rows):
 out=[];S=span(out)
 for r in rows:
  if tuple(r)not in S:out.append(np.array(r,dtype=np.int64));S=span(out)
 return np.array(out,dtype=np.int64)
def polarization(M,halfq=True):
 # M is the bilinear matrix over integers, quadratic is xMx/2 mod2.
 q=np.einsum('ij,jk,ik->i',V,M,V)//2%2
 t=[]
 for v,qq in zip(V,q):
  if qq or tuple(v)in span(t)or any(int(v@M@r)%2 for r in t):continue
  t.append(v.copy())
  if len(t)==4:break
 assert len(t)==4
 return np.array(t)
def complement(M,T):
 q=np.einsum('ij,jk,ik->i',V,M,V)//2%2;F=[]
 for j in range(4):
  target=np.eye(4,dtype=int)[j]
  candidates=[v for v,qq in zip(V,q)if qq==0 and np.array_equal(v@M@T.T%2,target)and all(int(v@M@r)%2==0 for r in F)]
  assert candidates;F.append(candidates[0].copy())
 F=np.array(F);assert len(span(list(T)+list(F)))==256
 assert not np.any(F@M@F.T%2)
 return F
def hnf(rows,n):
 H=fmpz_mat(np.array(rows,dtype=object).tolist()).hnf();a=np.array(H.tolist(),dtype=object);assert not np.any(a[n:]);return np.array(a[:n],dtype=np.int64)
def inverse_integral(A,B):
 Q=fmpq_mat(np.array(A,dtype=object).tolist())*fmpq_mat(np.array(B,dtype=object).tolist()).inv()
 assert all(x.denominator==1 for row in Q.tolist()for x in row)
 return np.array([[int(x)for x in row]for row in Q.tolist()],dtype=np.int64)
T1=np.array(alignment['polarization_basis'][:4],dtype=np.int64);T2=np.array(alignment['polarization_basis'][4:],dtype=np.int64)
assert len(span(list(T1)+list(T2)))==256 and np.array_equal(T1@G@T2.T%2,np.eye(4,dtype=int))
assert not np.any(T1@G@T1.T%2)and not np.any(T2@G@T2.T%2)
assert all(int(t@G@t)%4==0 for t in np.r_[T1,T2])
assert not np.any(np.einsum('ij,jk,ik->i',np.array(list(span(T1))),G,np.array(list(span(T1))))%4)
F=hnf(np.r_[2*np.eye(8,dtype=int),T1],8);F2=hnf(np.r_[2*np.eye(8,dtype=int),T2],8);G1=F@G@F.T
assert int(fmpz_mat(F.tolist()).det())==16 and not np.any(G1%2)and not np.any(np.diag(G1)%4)
assert int(fmpz_mat((G1//2).tolist()).det())==1
# Exact Lemma4.4 binary code.
C=np.array([[1,0,0,0,0,1,1,1,0,0],[0,1,0,0,0,0,1,1,1,0],[0,0,1,0,0,0,0,1,1,1],[0,0,0,1,0,1,0,0,1,1],[0,0,0,0,1,1,1,0,0,1]],dtype=np.int64)
Cp=np.c_[C[:,5:].T,np.eye(5,dtype=np.int64)];assert not np.any(C@Cp.T%2)
U=np.array(list(itertools.product(range(2),repeat=5)));CW=U@C%2;DW=U@Cp%2
for W in [CW,DW]:assert sorted(W.sum(1).tolist())==[0]+[4]*15+[6]*15+[10]
assert set(map(tuple,CW))&set(map(tuple,DW))=={tuple([0]*10),tuple([1]*10)}
rows=[2*r for r in np.eye(80,dtype=np.int64)]
for T,code in [(T1,C),(T2,Cp)]:
 for t in T:
  for c in code:rows.append(np.concatenate([int(a)*t for a in c]))
LC=hnf(rows,80);assert int(fmpz_mat(LC.tolist()).det())==2**40
sums=LC.reshape(80,10,8).sum(1);assert not np.any(sums%2)
constraint=(sums//2)@G@T1.T%2
# Find binary nullspace of constraint.T by independent Gaussian elimination.
M=constraint.T.copy();piv=[];r=0
for j in range(80):
 ids=np.flatnonzero(M[r:,j])
 if not len(ids):continue
 k=r+ids[0];M[[r,k]]=M[[k,r]]
 for i in range(4):
  if i!=r and M[i,j]:M[i]^=M[r]
 piv.append(j);r+=1
 if r==4:break
assert r==4
null=[]
for j in range(80):
 if j in piv:continue
 v=np.zeros(80,dtype=np.int64);v[j]=1
 for i,k in enumerate(piv):v[k]=M[i,j]
 assert not np.any(v@constraint%2);null.append(v)
K=hnf(np.r_[2*np.eye(80,dtype=int),null],80);Lp=K@LC
assert int(fmpz_mat(K.tolist()).det())==16
# Work in L1/2L1 with bilinear G1/2 and q=xG1x/4.
old=inverse_integral(2*np.eye(8,dtype=int),F);Aq=independent(old%2);assert len(Aq)==4
assert not np.any(Aq@(G1//2)@Aq.T%2)and all(int(v@G1@v)%8==0 for v in Aq)
Bq=complement(G1//2,Aq);Z=Bq@F
assert all(int(v@G@v)%8==0 for v in Z)
Ldouble2=hnf(np.r_[2*LC,np.array([np.tile(z,10)for z in F])],80)
assert int(fmpz_mat(Ldouble2.tolist()).det())==2**116
psi=Ldouble2.reshape(80,10,8).sum(1);assert not np.any(psi%2);psi=inverse_integral(psi//2,F)
# Twice root-coordinate basis for Lambda = L' + half-diagonal Z.
extra=np.array([np.tile(z,10)for z in Z]);Lambda2=hnf(np.r_[2*Lp,extra],80)
block=np.kron(np.eye(10,dtype=np.int64),G)
ldnorm=np.einsum('ij,jk,ik->i',Ldouble2,block,Ldouble2);assert not np.any(ldnorm%8)
assert np.array_equal(ldnorm//8%2,np.einsum('ij,jk,ik->i',psi,G1,psi)//4%2)
assert np.array_equal((Ldouble2@block@Ldouble2.T)//4%2,psi@(G1//2)@psi.T%2)
GG=Lambda2@block@Lambda2.T
assert not np.any(GG%8);GG//=8
assert not np.any(np.diag(GG)%2)and int(fmpz_mat(GG.tolist()).det())==1
assert int(fmpz_mat(Lambda2.tolist()).det())==2**120
assert np.array_equal(GG,GG.T)
att=np.zeros(80,dtype=np.int64);att[0]=2;att[8]=-2;attcoef=inverse_integral([2*att],Lambda2)[0]
assert int(attcoef@GG@attcoef)==8
# All definitions checked by membership in exact rational bases.
assert np.array_equal(inverse_integral(2*Lp,Lambda2)@Lambda2,2*Lp)
assert np.array_equal(inverse_integral(extra,Lambda2)@Lambda2,extra)
np.savez_compressed(D/'lattice80.npz',gram=GG,twice_root_coordinate_basis=Lambda2,e8_twice_euclidean_basis=E,e8_gram=G,LC_root_basis=LC,Lprime_root_basis=Lp,Ldoubleprime_twice_root_basis=Ldouble2,sum_map_in_L1_basis=psi,L1_root_basis=F,L2_root_basis=F2,T1=T1,T2=T2,binary_code=C,binary_dual=Cp,old_quotient_basis=Aq,new_quotient_basis=Bq,diagonal_lifts=Z,minimum8_coefficients=attcoef)
obj={'status':'PASS finite construction; minimum8 inherits primary Theorem4.7 subject to L1/sqrt2=E8 identification','dimension':80,'gram_determinant':'1','gram_even':True,'coordinate_convention':'Rows twice_root_coordinate_basis/2 are coefficients in10 copies of E8 root basis; metric is blockdiag(e8_gram)/2. Physical Euclidean rows equal Lambda2 blockdiag(E)/(4 sqrt2).','binary_weight_distribution':{'0':1,'4':15,'6':15,'10':1},'LC_index_in_E8power':str(2**40),'Lprime_index_in_LC':16,'Lambda_index_over_Lprime':16,'L1_determinant':int(fmpz_mat(G1.tolist()).det()),'L1_half_gram':(G1//2).tolist(),'T1':T1.tolist(),'T2':T2.tolist(),'old_quotient_basis':Aq.tolist(),'new_quotient_basis':Bq.tolist(),'quotient_pairing':(G1//2%2).tolist(),'quotient_quadratic_values':(np.einsum('ij,jk,ik->i',V,G1,V)//4%2).tolist(),'minimum8_coefficients':attcoef.tolist()}
(D/'construction.json').write_text(json.dumps(obj,indent=2)+'\n');print('PASS determinant1 even rank80 and minimum8 witness',flush=True);print('L1halfGram',G1//2,flush=True)
