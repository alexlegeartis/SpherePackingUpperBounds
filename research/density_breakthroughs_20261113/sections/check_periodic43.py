"""Self-contained replay of the compact numerical43D contact-weight witness.

Only periodic43_fista.npz is required, not the external clone or ambient pool.
Exact contact/saturation assertions are separate from numerical residual/rank.
"""
from pathlib import Path
import json,time
import numpy as np
import sympy as s
from sympy.matrices.normalforms import smith_normal_form
P=Path(__file__).resolve().parent;start=time.time()
d=np.load(P/'periodic43_fista.npz')
G=d['G'].astype(np.int64);T=d['T'].astype(np.int64);V=d['V'].astype(np.int64)
pairs=d['pairs'].astype(np.int64);E=d['E'];weights=d['weights']
assert G.shape==(48,48) and T.shape==(5,48) and V.shape==(12000,48)
assert pairs.shape==(len(V),2) and np.all((0<=pairs)&(pairs<6))
assert weights.shape==(len(V),) and np.all(weights>=0)
K=T@G@T.T;assert int(s.Matrix(K).det())==3072
sm=smith_normal_form(s.Matrix(T),domain=s.ZZ);assert all(abs(sm[i,i])==1 for i in range(5))
adj=np.array(s.Matrix(K).inv()*3072,dtype=np.int64)
assert np.all(np.einsum('ij,ij->i',V@G,V)==6)
S=np.vstack([np.zeros(5,dtype=np.int64),np.eye(5,dtype=np.int64)])
F=S[pairs[:,0]]-S[pairs[:,1]]
assert np.array_equal(V@G@T.T,F)
assert np.all(np.einsum('ij,ij->i',F@adj,F)==896)
L=np.linalg.cholesky(G)
assert E.shape==(48,43) and np.linalg.norm(E.T@E-np.eye(43))<1e-10
assert np.linalg.norm((T@L)@E)<1e-10
W=V@L@E
assert np.max(abs(np.einsum('ij,ij->i',W,W)-137/24))<1e-8
n=43;ii,jj=np.triu_indices(n);q=len(ii)
A=np.zeros((len(V),q+5*n));A[:,:q]=W[:,ii]*W[:,jj]*(1+(ii!=jj))
for i in range(1,6):A[:,q+(i-1)*n:q+i*n]=2*W*((pairs[:,0]==i).astype(int)-(pairs[:,1]==i).astype(int))[:,None]
goal=np.zeros(A.shape[1]);goal[:q]=(ii==jj)
residual=float(np.linalg.norm(A.T@weights-goal));assert residual<1e-8
eig=np.linalg.eigvalsh(A.T@A);rank=int(np.sum(eig>eig[-1]*1e-10));assert rank==1161
report={'exact_checks':'all12000 ambient norms6, oriented normal differences, squared normal norms7/24, normal determinant3072 and saturation', 'numerical_rank':rank,'variables':A.shape[1],'positive_weights':int(np.sum(weights>1e-12)),'nonnegative_weights':bool(np.all(weights>=0)),'trace_representation_residual':residual,'elapsed':time.time()-start,'scope':'numerical first-order obstruction candidate; not an exact eutaxy certificate or a nonlinear local-optimality theorem'}
(P/'periodic43_witness_checked.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
