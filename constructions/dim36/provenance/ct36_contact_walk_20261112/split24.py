import numpy as np,sympy as sp,json
from pathlib import Path
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parent;S=P.parent/'ct36_metric_screen_20261112';H=np.load(S/'data.npz')['H'];V=np.load(S/'min6.npy');A=V.T@V@H;M=sp.Matrix(A@(A-192*np.eye(36,dtype=np.int64)));D,U,W=smith_normal_decomp(DomainMatrix.from_Matrix(M).convert_to(sp.ZZ));B=hermite_normal_form(W.to_Matrix()[:,12:]).T;assert B.shape==(24,36);assert M*B.T==sp.zeros(36,24)
T=sp.Matrix(np.load(P.parent/'ct36_symmetry_break_20261112/reduced.npz')['T']);C=B*T.inv();assert all(x.q==1 for x in C);C=np.array(C,dtype=np.int64);B12=np.load(P/'block_basis.npy')[12:24];full=np.vstack([B12,C]);index=abs(sp.Matrix(full).det());np.save(P/'split_basis.npy',full);print('index',index,flush=True);json.dump({'spectator_rank':12,'coupled_rank':24,'index':str(index)},open(P/'split24.json','w'))
