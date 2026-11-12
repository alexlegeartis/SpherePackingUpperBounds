import numpy as np,sympy as s,json
from pathlib import Path
P=Path('research/ct36_metric_screen_20261112');x=np.load('research/ct36_contact_walk_20261112/flex.npz');T=np.load('research/ct36_symmetry_break_20261112/reduced.npz')['T'];I=np.array(s.Matrix(T.tolist()).inv()).astype(np.int64)
for file,key in [('min6.npy',None),('newshort.npy',None),('kernel_sat_basis.npy',None),('types.npz','representatives')]:
 V=np.load(P/file);V=V[key] if key else V;V=V.reshape(-1,36)@I
 Q=np.einsum('bi,ij,bj->b',V,x['G4'],V);R=np.einsum('bi,ij,bj->b',V,x['D'],V)
 print(file,'minG',Q.min()/4,'minminus',min(Q/4-R/5),'minplus',min(Q/4+R/5))
 if np.any((Q==32)&(R!=0)): print('BAD CONTACT',V[np.where((Q==32)&(R!=0))[0][0]])
