"""Sampled radial moment-cone QP. No value is a certificate by itself."""
from pathlib import Path
from fractions import Fraction as Q
import json,math,numpy as np,cvxpy as cp

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'verifications/dim4/refinement_energy20261115/vendor/D4/data'
k=json.loads((DATA/'v7/kernel.json').read_text());o=json.loads((DATA/'v7/origin.json').read_text())
rays=k['blocks'][0]['rays'];den=int(k['common_denominator'])
A=np.array([[float(Q(x)) for x in row['vector'][1:]] for row in rays])
u=np.array([float(Q(x)) for x in o['u']]);b=A@u
ww=np.sqrt([float(Q(row['weight'],den)) for row in rays])
A*=ww[:,None]*1000;b*=ww*1000
breaks=list(map(Q,k['breaks']));intervals=[tuple(map(Q,row)) for row in k['radial_basis_intervals']]
features=[];pieces=[];radii=[]
for i in range(7):
    left,right=intervals[i]
    for radius in np.linspace(float(breaks[i]),float(breaks[i+1]),65):
        tau=(radius-float(left))/float(right-left)
        phi=np.zeros(49)
        phi[7*i:7*i+7]=[math.comb(6,j)*tau**j*(1-tau)**(6-j) for j in range(7)]
        features.append(phi);pieces.append(i);radii.append(radius)
F=np.array(features).T;pieces=np.array(pieces)
profile=json.loads((DATA/'census_floors.json').read_text())['profile']
results=[]
for key in [(24,1),(24,2),(25,4),(25,5),(26,5)]:
    row=next(g for g in profile if (g['N'],g['first_radius'])==key)
    counts=row['minimizing_counts'];diffs=np.diff([0]+counts)
    weights=cp.Variable(len(pieces),nonneg=True)
    constraints=[cp.sum(weights[pieces==i])==diffs[i] for i in range(5)]
    # The elementary volume bound gives at most 255 nonorigin centres in
    # radius 3/2. The final piece is a subset of the inherited 151-point
    # annulus [33/28,7/4]. These are relaxed, valid upper bounds.
    constraints += [cp.sum(weights[pieces==5])<=255-counts[-1],cp.sum(weights[pieces==6])<=151]
    problem=cp.Problem(cp.Minimize(cp.sum_squares(A@F@weights+b)),constraints)
    try:problem.solve(solver='CLARABEL',tol_gap_abs=1e-8,tol_feas=1e-9,tol_gap_rel=1e-9,max_iter=250,max_threads=1)
    except cp.error.SolverError:problem.solve(solver='SCS',eps=1e-8,max_iters=100000)
    result=dict(N=key[0],j=key[1],counts=counts,status=problem.status,E0_candidate=problem.value/1e6,
                z=(F@weights.value).tolist(),active_radii=[(radii[i],float(x)) for i,x in enumerate(weights.value) if x>1e-4])
    results.append(result)
    print({kk:vv for kk,vv in result.items() if kk!='z'},flush=True)
(Path(__file__).resolve().parent/'radial_moment_qp_screen.json').write_text(json.dumps(results,indent=2)+'\n')
