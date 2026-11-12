import numpy as np,sympy as sp,json,time
from pathlib import Path
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parent;S=P.parent/'ct36_metric_screen_20261112';H=np.load(S/'data.npz')['H'];V=np.load(S/'min6.npy');A=V.T@V@H;Bs=[]
for l in [0,96,192]:
 M=sp.Matrix(A-l*np.eye(36,dtype=np.int64));S,U,W=smith_normal_decomp(DomainMatrix.from_Matrix(M).convert_to(sp.ZZ));B=hermite_normal_form(W.to_Matrix()[:,24:]).T;assert B.shape==(12,36);assert M*B.T==sp.zeros(36,12);Bs.append(B);print('lambda',l,'gramdet',str((B*sp.Matrix(H)*B.T).det()),flush=True)
B=sp.Matrix.vstack(*Bs);T=sp.Matrix(np.load(P.parent/'ct36_symmetry_break_20261112/reduced.npz')['T']);C=B*T.inv();assert all(x.q==1 for x in C);np.save(P/'block_basis.npy',np.array(C,dtype=np.int64));print('index',B.det(),flush=True);json.dump({'index':str(abs(B.det()))},open(P/'decompose.json','w'))
