#!/usr/bin/env python3
"""Reconstruct the one-dimensional 56D contact nullspace exactly."""
from pathlib import Path
import json,math,time
import numpy as np
import sympy as sp
from scipy.linalg import eigh

P=Path(__file__).resolve().parent
d=np.load(P/'L56a_orbit.npz');G=d['G'];V=d['V'];n=len(G)
rng=np.random.default_rng(560013);L=np.linalg.cholesky(G)
W=V[rng.choice(len(V),6000,False)]@L
i,j=np.triu_indices(n);A=W[:,i]*W[:,j]*np.where(i==j,1,2)
ev,vec=eigh(A.T@A,subset_by_index=[0,0]);Dp=np.zeros((n,n));Dp[i,j]=vec[:,0];Dp[j,i]=vec[:,0]
Df=L@Dp@L.T;Df/=np.max(abs(Df));reconstructed=None
for scale in range(1,10001):
    D=np.rint(Df*scale).astype(np.int64)
    if np.max(abs(D-Df*scale))>1e-6:continue
    if np.any(np.einsum('bi,ij,bj->b',V,D,V)):continue
    reconstructed=D;break
assert reconstructed is not None
D=reconstructed
content=math.gcd(*D.ravel().tolist());D//=content
if D[np.unravel_index(np.argmax(abs(D)),D.shape)]<0:D=-D
assert np.all(np.einsum('bi,ij,bj->b',V,D,V)==0)
Gs=sp.Matrix(G);Ds=sp.Matrix(D);operator=Gs.inv()*Ds
print('maxD',int(abs(D).max()),'trace',operator.trace(),'trace2',(operator*operator).trace(),flush=True)
# A quadratic minimal polynomial is expected from the real quadratic structure.
a=sp.simplify(sp.trace(operator)/n);b=sp.simplify(sp.trace(operator*operator)/n)
scalar_square=operator*operator==b*sp.eye(n)
np.savez_compressed(P/'L56a_flex.npz',G=G,D=D,V=V)
report={'dimension':n,'contact_vectors_up_to_sign':len(V),'reconstruction_scale':scale,'D_max_abs':int(abs(D).max()),'exact_all_contact_slopes_zero':True,'exact_operator_trace':str(operator.trace()),'exact_operator_square_scalar':bool(scalar_square),'operator_square_scalar':str(b),'determinant_if_scalar_square':f'29^28 * (1 - ({b}) * t^2)^28' if scalar_square else None,'scope':'exact contact null direction; other lattice vectors are not yet constrained'}
(P/'L56a_flex.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
np.savetxt(P/'L56a_D.txt',D,fmt='%d')
