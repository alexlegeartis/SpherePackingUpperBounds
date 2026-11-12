import numpy as np,json,sympy as sp
from fractions import Fraction as F
from pathlib import Path
P=Path(__file__).resolve().parent;S=P.parent/'ct36_symmetry_break_20261112';V=np.loadtxt(S/'vectors.txt',dtype=np.int64).reshape(-1,36);G=np.load(S/'reduced.npz')['N']/4;i,j=np.triu_indices(36);rng=np.random.default_rng(345);v=V[rng.choice(len(V),2400,False)];A=v[:,i]*v[:,j];_,sv,vh=np.linalg.svd(A.astype(float),full_matrices=False);u=vh[-1];u/=max(abs(u));fr=[F(float(x)).limit_denominator(10000)for x in u];den=int(sp.ilcm(*[x.denominator for x in fr]));z=np.array([int(x*den)for x in fr],dtype=np.int64);print('rank',sum(sv>1e-7),'den',den,'max',max(abs(z)),flush=True)
# upper coefficient form, offdiagonal coefficients are twice Gram entries
D=np.zeros((36,36),dtype=np.int64);D[i,j]=z;D[j,i]=z;D[np.diag_indices(36)]*=2
q=np.einsum('bi,ij,bj->b',V,D,V);print('constraints',len(V),'maxresid',max(abs(q)),flush=True);assert not q.any()
np.savez(P/'flex.npz',G4=np.load(S/'reduced.npz')['N'],D=D,V=V)
A=np.linalg.solve(G,D);ev=np.linalg.eigvals(A).real;print('trace',np.trace(A),'ev',ev.min(),ev.max(),flush=True)
json.dump({'sample_rank':int(sum(sv>1e-7)),'constraints_verified':len(V),'direction_denominator':den,'trace_float':float(np.trace(A)),'eigen_range':[float(ev.min()),float(ev.max())]},open(P/'flex.json','w'),indent=2)
