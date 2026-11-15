"""Numerical diagnostics for finite-count radial Gram projection floors."""
from pathlib import Path
from fractions import Fraction as Q
import json,numpy as np,cvxpy as cp

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'verifications/dim4/refinement_energy20261115/vendor/D4/data'
j=json.loads((DATA/'v7/kernel.json').read_text())
den=int(j['common_denominator'])
rays=j['blocks'][0]['rays']
A=np.array([[float(Q(x)) for x in row['vector'][1:]] for row in rays])
origin=json.loads((DATA/'v7/origin.json').read_text())
u=np.array([float(Q(x)) for x in origin['u']])
# The current origin is the separately saturated vector u. The first entry
# of the historical 50-coordinate ray is not the current origin feature.
b=A@u
weights=np.sqrt([float(Q(row['weight'],den)) for row in rays])
A*=weights[:,None]*1000;b*=weights*1000
profile=json.loads((DATA/'census_floors.json').read_text())['profile']
print('breaks:',j['breaks'],flush=True)
results=[]
for key in [(24,1),(24,2),(25,4),(25,5),(26,5)]:
    row=next(g for g in profile if (g['N'],g['first_radius'])==key)
    counts=row['minimizing_counts']
    diffs=np.diff([0]+counts)
    # Finite caps below are hypothetical diagnostics, not proved packing
    # constraints. Only the uncapped case is an unconditional relaxation.
    for cap in [None,100,70,60]:
        z=cp.Variable(49,nonneg=True)
        constraints=[cp.sum(z[7*k:7*k+7])==diffs[k] for k in range(5)]
        if cap:
            constraints += [cp.sum(z[35:42])<=cap-counts[-1],cp.sum(z[42:])<=151]
        problem=cp.Problem(cp.Minimize(cp.sum_squares(A@z+b)),constraints)
        try:
            problem.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_feas=1e-9,tol_gap_rel=1e-9,max_iter=250,max_threads=1)
        except cp.error.SolverError:
            problem.solve(solver='SCS',eps=1e-8,max_iters=50000)
        if z.value is None:continue
        result=dict(N=key[0],j=key[1],counts=counts,cap=cap,additional_cap_proved=False,certified=False,status=problem.status,
                    E0_candidate=problem.value/1e6,outer_totals=[float(sum(z.value[35:42])),float(sum(z.value[42:]))],
                    z=z.value.tolist())
        results.append(result)
        print({k:v for k,v in result.items() if k!='z'},flush=True)
(Path(__file__).resolve().parent/'census_qp_screen.json').write_text(json.dumps(results,indent=2)+'\n')
