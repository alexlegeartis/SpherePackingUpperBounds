import numpy as np,sympy as s,json,sys
from pathlib import Path
p=Path(__file__).parent;f=np.load(p/'bkz36.npz');quarter=len(sys.argv)>1;U=s.Matrix(f['U']);base=s.Matrix(np.load(p.parent/'ct36_contact_walk_20261112'/('candidate_quarter.npy' if quarter else 'candidate_1.npy')));Q=U*base*U.T
assert abs(U.det())==1 and Q==U*base*U.T
assert all(int(Q[i,i])%2==0 for i in range(36))
L,D=Q.LDLdecomposition(hermitian=False);assert L*D*L.T==Q and all(D[i,i]>0 for i in range(36))
S=2**40
with open(p/('ldl_quarter.txt' if quarter else 'ldl_bounds.txt'),'w')as w:
 for i in range(36):
  for x in [D[i,i]]+[L[j,i]for j in range(36)]:
   lo=int(s.floor(x*S));hi=int(s.ceiling(x*S));assert abs(lo)<2**53 and abs(hi)<2**53;w.write(f'{lo} {hi}\n')
json.dump({'scale':S,'dimension':36,'bound':159,'even_integral':True,'unimodular':True,'exact_ldl_identity':True},open(p/('setup_quarter.json' if quarter else 'setup.json'),'w'),indent=2)
print('Exact LDL interval inputs and unimodular basis checked.')
