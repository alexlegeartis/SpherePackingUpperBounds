"""Audit compact exact contacts and numerical periodic eutaxy weights.

This does not certify the minimum of the 48D mother lattice or turn floating
eutaxy weights into an exact local-optimality theorem.
"""
import json
import pathlib
import sys
import re
import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from scipy.linalg import null_space

p=pathlib.Path(sys.argv[1]) if len(sys.argv)>1 else pathlib.Path(__file__).resolve().parent/'dim48mod3'
compact=p/'periodic44_eutaxy_sample.npz'
if not compact.exists():
    state=np.load(p/'periodic44_adaptive_state.npz')
    e=np.load(p/'antipode44_embedding.npz')
    V=state['V'][state['selected']]
    assert np.max(np.abs(V))<32768
    data=dict(V=V.astype(np.int16),weights=state['weights'],P=e['P'],N=e['N'],w=e['w'])
    if 'source_name' in e:data['source_name']=e['source_name']
    np.savez_compressed(compact,**data)
f=np.load(compact);V=f['V'].astype(np.int64);P,N,w=f['P'],f['N'],f['w'];weights=f['weights']
assert sp.Matrix(P.tolist()).det()==1
assert np.all(P==P.T) and np.all(P.diagonal()%2==0)
G=N@P@N.T;assert sp.Matrix(G.tolist()).det()==605
smith=smith_normal_form(sp.Matrix(N.tolist()),domain=sp.ZZ)
invariants=[abs(int(smith[i,i])) for i in range(4)]
assert invariants==[1,1,1,1]
BN=np.array((605*sp.Matrix(G.tolist()).inv()).tolist(),np.int64)
F=V@w.T
assert np.array_equal(w,N@P)
assert np.all(np.einsum('vi,ij,vj->v',V,P,V)==6)
assert np.all(np.einsum('vi,ij,vj->v',F,BN,F)==176)
S=np.vstack([np.zeros((1,4),int),np.eye(4,dtype=int)])
pairs={tuple(S[i]-S[j]):(i,j) for i in range(5) for j in range(5) if i!=j}
idx=np.array([pairs[tuple(f)] for f in F])
if 'pairs' in f:assert np.array_equal(idx,f['pairs'])
C=np.linalg.cholesky(P);Z=null_space(N@C);Y=V@C@Z
ii,jj=np.triu_indices(44);nmetric=len(ii)
A=np.zeros((len(Y),nmetric+176));A[:,:nmetric]=Y[:,ii]*Y[:,jj]*(1+(ii!=jj))
for i in range(1,5):
    signs=(idx[:,0]==i).astype(float)-(idx[:,1]==i)
    A[:,nmetric+44*(i-1):nmetric+44*i]=2*Y*signs[:,None]
goal=np.zeros(A.shape[1]);goal[:nmetric][ii==jj]=1
assert np.all(weights>=0)
positive=weights>1e-8
eig=np.linalg.eigvalsh(A[positive].T@A[positive])
# These are explicit numerical acceptance tolerances, not exact proofs.
residual_2=float(np.linalg.norm(A.T@weights-goal))
residual_inf=float(np.max(np.abs(A.T@weights-goal)))
support_rank=int((eig>eig[-1]*1e-10).sum())
assert residual_2<1e-7 and residual_inf<1e-7
assert support_rank==1166
minimum_status='minimum6 not certified here; exact norm6 vectors exhibited'
if 'source_name' in f and str(f['source_name'])=='CQ48a':
    source=(pathlib.Path(__file__).resolve().parent/'sources'/'CQ48a.html').read_text()
    section=re.search(r'<a\s+NAME="GRAM[^\"]*">.*?</a>(.*?)(?=<p><li>|<a\s+NAME=|<hr>|$)',source,re.S|re.I).group(1)
    values=list(map(int,re.findall(r'-?\d+',re.sub(r'<[^>]+>',' ',section))))
    if values[:2]==[48,48] and len(values)==2306:
        archived=np.array(values[2:]).reshape(48,48)
    elif values[:2]==[48,0] and len(values)==1178:
        archived=np.zeros((48,48),np.int64);a,b=np.tril_indices(48)
        archived[a,b]=values[2:];archived[b,a]=values[2:]
    else:raise ValueError('Unexpected archived catalogue Gram format')
    assert np.array_equal(P,archived)
    minimum_status='published minimum6 inherited from catalogue CQ48a (P48n); exact Gram identity checked'
report={'contacts':len(V),'variables':A.shape[1],
        'mother_gram_determinant_exact':1,'mother_even_integral_exact':True,
        'mother_minimum_status':minimum_status,
        'normal_gram_determinant_exact':605,'contact_projected_squared_norm_exact':'314/55',
        'normal_smith_invariants_exact':invariants,'normal_lattice_primitive_exact':True,
        'weight_minimum':float(weights.min()),'positive_weights_above_1e-8':int(positive.sum()),
        'smallest_positive_support_gram_eigenvalue':float(eig[0]),
        'positive_support_rank_tolerance_1e-10':support_rank,
        'eutaxy_residual_2norm':residual_2,
        'eutaxy_residual_infinity':residual_inf,
        'numerical_residual_acceptance_tolerance':1e-7,
        'scope':'exact integer contact checks; numerical rank and eutaxy only; no new density record'}
(p/'periodic44_eutaxy_audit.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2),flush=True)
