"""Numerical tangent screen for actual five-coset distances."""
import json
import pathlib
import numpy as np
from scipy.optimize import linprog

p = pathlib.Path(__file__).resolve().parent / 'dim48mod3'
e = np.load(p / 'antipode44_embedding.npz')
P, w = e['P'], e['w']
f = np.load(p / 'flexes.npz')
D = f['Ds'] * (3. / (4 * 206))
V = np.load(p / 'periodic44_initial.npz')['V']
for source in sorted(p.glob('periodic44_shell*.npy')):
    V = np.vstack([V, np.load(source)])
V = np.unique(V, axis=0)
F = V @ w.T
Pi = np.linalg.inv(P)
M = w @ Pi @ w.T
B = np.linalg.inv(M)
normal = F @ B @ w @ Pi
A = (np.einsum('vi,aij,vj->va', V, D, V)
     - np.einsum('vi,aij,vj->va', normal, D, normal))
m = 314./55
grad = -0.5 * np.array([np.trace(Pi @ d) - np.trace(B @ w @ Pi @ d @ Pi @ w.T) for d in D])
grad = np.append(grad,22/m)
constraints = np.column_stack([-A, np.ones(len(A))])
res = linprog(-grad, A_ub=constraints, b_ub=np.zeros(len(A)),
              bounds=[(-1,1)]*7,method='highs')
report = {'contacts':len(V), 'unique_derivative_constraints':len(np.unique(np.rint(A*1e7).astype(np.int64),axis=0)),
          'success':bool(res.success),'objective_improvement':float(-res.fun),
          'direction':res.x.tolist(),'gradient':grad.tolist(),
          'scope':'floating tangent-cone screen on sampled exact contacts'}
print(json.dumps(report,indent=2),flush=True)
(p / 'periodic44_tangent_report.json').write_text(json.dumps(report,indent=2))
np.savez(p / 'periodic44_tangent.npz', V=V, A=A, gradient=grad, direction=res.x)
