"""Numerical full periodic tangent screen for the five-coset 44D packing.

Each input is an integer vector of exact P-norm 6 with one of five active
normal differences.  Variables are a symmetric 44D metric and four
independent 44D coset displacements.  This is a discovery calculation;
the sampled contact set is not asserted complete.
"""
import json
import pathlib
import sys
import time
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import linprog

p = pathlib.Path(__file__).resolve().parent / 'dim48mod3'
e = np.load(p / 'antipode44_embedding.npz')
P, N, w = e['P'], e['N'], e['w']
V = np.load(p / 'periodic44_initial.npz')['V']
for source in sorted(p.glob('periodic44_shell*.npy')):
    V = np.vstack([V, np.load(source)])
V = np.unique(V,axis=0)
rng = np.random.default_rng(42)
count = min(len(V),int(sys.argv[1]) if len(sys.argv)>1 else 6000)
V = V[rng.choice(len(V),count,replace=False)]
assert np.all(np.einsum('ij,jk,ik->i',V,P,V)==6)
F = V @ w.T
S = np.vstack([np.zeros((1,4),np.int64),np.eye(4,dtype=np.int64)])
pairs = {tuple(S[i]-S[j]):(i,j) for i in range(5) for j in range(5) if i!=j}
C = np.linalg.cholesky(P)
Z = null_space(N @ C)
Y = V @ C @ Z
m = 314./55
assert np.max(np.abs(np.einsum('ij,ij->i',Y,Y)-m))<1e-9
ii,jj = np.triu_indices(44)
M = Y[:,ii]*Y[:,jj]*(1+(ii!=jj))
T = np.zeros((len(Y),176))
for k,(y,f) in enumerate(zip(Y,F)):
    i,j = pairs[tuple(f)]
    if i:T[k,44*(i-1):44*i] += 2*y
    if j:T[k,44*(j-1):44*j] -= 2*y
A = np.column_stack([M,T])
goal = np.zeros(A.shape[1]); goal[:len(ii)][ii==jj]=1
eig,U = np.linalg.eigh(A.T @ A)
keep=eig>eig[-1]*1e-10
t = U[:,keep]@((U[:,keep].T@goal)/eig[keep])
weights=A@t
report={'contacts':len(V),'variables':A.shape[1],
        'numerical_rank':int(keep.sum()),'gram_smallest_eigenvalues':eig[:12].tolist(),
        'least_norm_weights_min':float(weights.min()),
        'least_norm_weights_negative':int((weights<0).sum()),
        'least_norm_goal_residual':float(np.linalg.norm(A.T@weights-goal)),
        'scope':'sampled numerical tangent-cone calculation; no new packing claim'}
print(json.dumps(report,indent=2),flush=True)
if len(sys.argv)>2 and sys.argv[2]=='fista':
    X=(A@U[:,keep])/np.sqrt(eig[keep])
    target=(U[:,keep].T@goal)/np.sqrt(eig[keep])
    best=np.maximum(weights,0); y=best.copy(); t=1.; began=time.time()
    for iteration in range(20001):
        residual=X.T@y-target
        nxt=np.maximum(0,y-X@residual)
        nt=(1+np.sqrt(1+4*t*t))/2
        y=nxt+(t-1)/nt*(nxt-best)
        best=nxt; t=nt
        if iteration%100==0:
            err=float(np.linalg.norm(X.T@best-target))
            if iteration%1000==0:
                print('FISTA',iteration,'residual',err,'seconds',time.time()-began,flush=True)
            if err<1e-10:break
    report['fista_iterations']=iteration
    report['fista_transformed_residual']=err
    report['fista_original_residual']=float(np.linalg.norm(A.T@best-goal))
    report['fista_nonzero_weights']=int((best>1e-12).sum())
    weights=best
    print('FISTA final',json.dumps({k:v for k,v in report.items() if k.startswith('fista')}),flush=True)
elif weights.min() < -1e-10 or np.linalg.norm(A.T@weights-goal)>1e-7:
    res=linprog(goal,A_ub=-A,b_ub=np.zeros(len(A)),bounds=[(-1,1)]*len(goal),method='highs')
    report['lp_success']=bool(res.success)
    report['minimum_trace']=float(res.fun) if res.success else None
    report['lp_message']=res.message
    if res.success:
        np.savez(p/'periodic44_full_direction.npz',direction=res.x,Z=Z,C=C,ii=ii,jj=jj,V=V,Y=Y,A=A)
    print('LP',json.dumps({k:v for k,v in report.items() if k.startswith('lp') or k=='minimum_trace'}),flush=True)
else:
    report['positive_eutaxy_on_sample']=True
np.savez(p/'periodic44_full_cone.npz',V=V,Y=Y,weights=weights,goal=goal,rank=int(keep.sum()))
(p/'periodic44_full_cone_report.json').write_text(json.dumps(report,indent=2))
