"""Exact primitive period basis and five rational cosets for the 44D family."""
import json
import pathlib
import numpy as np
import sympy as sp
from fpylll import IntegerMatrix, BKZ

p=pathlib.Path(__file__).resolve().parent/'dim48mod3'
e=np.load(p/'antipode44_embedding.npz');P=e['P'];w=e['w']
A=sp.Matrix(w.tolist());U=sp.eye(48)
for k in range(4):
    nz=next(j for j in range(k,48) if A[k,j])
    A.col_swap(k,nz);U.col_swap(k,nz)
    for j in range(k+1,48):
        a,b=A[k,k],A[k,j]
        if not b:continue
        x,y,g=sp.gcdex(a,b)
        ck,cj=A[:,k],A[:,j];uk,uj=U[:,k],U[:,j]
        A[:,k]=x*ck+y*cj;A[:,j]=(-b//g)*ck+(a//g)*cj
        U[:,k]=x*uk+y*uj;U[:,j]=(-b//g)*uk+(a//g)*uj
    if A[k,k]==-1:A[:,k]=-A[:,k];U[:,k]=-U[:,k]
    assert A[k,k]==1
    for j in range(k):
        q=A[k,j];A[:,j]-=q*A[:,k];U[:,j]-=q*U[:,k]
assert A==sp.eye(4).row_join(sp.zeros(4,44))
T=np.array(U.T.tolist(),dtype=np.int64)
period=T[4:]
G=period@P@period.T
R=np.rint(np.linalg.cholesky(G)*1e8).astype(np.int64)
L=IntegerMatrix.from_matrix(R.tolist());V=IntegerMatrix.identity(44)
BKZ.reduction(L,BKZ.Param(block_size=25,max_loops=2),U=V,float_type='long double')
V=np.array([[int(V[i,j]) for j in range(44)] for i in range(44)],np.int64)
T[4:]=V@period;period=T[4:];G=period@P@period.T
cross=T[:4]@P@period.T
positions=sp.Matrix(cross.tolist())*sp.Matrix(G.tolist()).inv()
den=int(sp.ilcm(*[v.q for v in positions]))
Rnum=np.array((positions*den).tolist(),np.int64)
assert sp.Matrix(G.tolist()).det()==605
assert abs(sp.Matrix(T.tolist()).det())==1
assert np.array_equal(w@T.T,np.column_stack([np.eye(4,dtype=np.int64),np.zeros((4,44),np.int64)]))
np.savez(p/'periodic44_period.npz',T=T,period=period,G=G,position_numerators=Rnum,position_denominator=den)
(p/'periodic44_period.json').write_text(json.dumps({'period_gram_det':605,'positions_denominator':den,'period_gram':G.tolist(),'position_numerators':Rnum.tolist(),'scope':'exact algebraic period extraction; minimum not certified here'},indent=2))
print('det605 coset denominator',den,'diag',G.diagonal(),flush=True)
