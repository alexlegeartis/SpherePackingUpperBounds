#!/usr/bin/env python3
"""Exact rational error audit for coset_enum_logged.cpp.

This checks coefficient errors and a uniform pruning safety margin for every
integer vector of exact norm <=60. It does not substitute a completed traversal.
See NE64_CVP_CERTIFICATE.md for the proof connecting the bounds to the loop.
"""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import hashlib,json,numpy as np
P=Path(__file__).resolve().parent
N=64
D=np.load(P/'ne64_primal_best.npz')
G=D['G'].astype(object);H=D['H'].astype(object);T=D['T'].astype(object);Ti=D['Ti'].astype(object);Q=D['R'].astype(object)
I=np.eye(N,dtype=object)
assert np.array_equal(T@Ti,I) and np.array_equal(Ti@T,I)
assert np.array_equal(np.loadtxt(P/'ne64_primal_best_Q.txt',dtype=np.int64).astype(object),Q)
assert np.array_equal(T@G@T.T,Q)
assert np.array_equal(G@H,3*I)
Qinv3=Ti.T@H@Ti
assert np.array_equal(Q@Qinv3,3*I)

def sqrt_ceil(x):
    x=F(x)
    r=isqrt(x.numerator//x.denominator)
    return r if r*r*x.denominator==x.numerator else r+1

def enc(x):
    x=F(x)
    return {'rational':str(x),'decimal':float(x)}

# Rational LDL: Q=L diag(d) L^T.
L=[[F(int(i==j)) for j in range(N)] for i in range(N)]
d=[]
for i in range(N):
    for j in range(i):
        L[i][j]=(F(int(Q[i,j]))-sum((L[i][k]*d[k]*L[j][k] for k in range(j)),F(0)))/d[j]
    d.append(F(int(Q[i,i]))-sum((L[i][k]**2*d[k] for k in range(i)),F(0)))
    assert d[-1]>0
U=[[L[j][i] if j>=i else F(0) for j in range(N)] for i in range(N)]
hexpath=P/'ne64_certified_run_ldl.hex'
lines=hexpath.read_text().splitlines()
assert len(lines)==N
logged=[[F.from_float(float.fromhex(s)) for s in line.split()] for line in lines]
assert all(len(row)==N+1 for row in logged)
df=[row[0] for row in logged]
uf=[row[1:] for row in logged]
assert all(x>0 for x in df)
# The stored U diagonal is unused and zero. Only upper triangular entries matter.
assert all(uf[i][j]==0 for i in range(N) for j in range(i+1))
M=[sqrt_ceil(F(60*int(Qinv3[i,i]),3)) for i in range(N)]
Z=[sqrt_ceil(F(60)/x) for x in d]
u=F(1,2**53)
gamma=lambda k:k*u/(1-k*u)
# A center summand follows fewer than70 rounded operations. The gamma70
# bound is deliberately more conservative than the actual4-accumulator loop.
s=[sum((abs(uf[i][j])*M[j] for j in range(i+1,N)),F(0)) for i in range(N)]
eta=[sum((abs(uf[i][j]-U[i][j])*M[j] for j in range(i+1,N)),F(0))+gamma(70)*s[i] for i in range(N)]
Chat=[(1+gamma(70))*x for x in s]
# Eps is3 times an absolute error bound on fl((center-residue)/3).
eps=[gamma(2)*(x+2) for x in Chat]
B=[F(3,2)+x for x in eps]
assert max(eps)<F(1,100)
assert max((1+u)*b for b in B)<F(8,5)
component=[];base_error=[]
for i in range(N):
    # A step v is formed from fl(base-center) plus an exact integer3*step.
    kappa=u*B[i]+u*(Z[i]+eta[i]+u*B[i])
    error=eta[i]+kappa
    term=abs(df[i]-d[i])*Z[i]**2 + df[i]*(2*Z[i]*error+error**2) + gamma(2)*df[i]*(Z[i]+error)**2
    component.append(term)
    # Rounding the quotient can choose the neighboring integer only at a tie;
    # its squared-distance penalty is at most6*eps+4*eps^2.
    base=abs(df[i]-d[i])*Z[i]**2+df[i]*(2*Z[i]*eta[i]+eta[i]**2)+df[i]*(6*eps[i]+4*eps[i]**2)+((1+u)**4-1)*df[i]*B[i]**2
    base_error.append(base)
component_sum=sum(component,F(0))
partial_error=component_sum+gamma(64)*(60+component_sum)
# fl(R-partial) has absolute error <=u*61 along a true norm<=60 path.
remaining_error=61*u
radius=F.from_float(60.000001)
margin=radius-60
# At a true-coordinate branch, component errors are disjoint; using all
# component errors again is conservative and makes the audit simpler.
branch_bound=partial_error+max(component)+remaining_error
base_bound=partial_error+max(base_error)+remaining_error
assert partial_error<F(1,100)
assert max(branch_bound,base_bound)<margin

# Coarse arithmetic-domain caps for ALL visited nodes, independent of the
# desired-vector bounds M. These establish exact integer conversions and
# exclude signed64/__int128 overflow. The extra unit dominates rounding.
visited=[0]*N
for i in range(N-1,-1,-1):
    visited[i]=sqrt_ceil(F(61)/df[i])+sqrt_ceil(sum((abs(uf[i][j])*visited[j] for j in range(i+1,N)),F(0))**2)+1
visited_sumabs=[sum((abs(uf[i][j])*visited[j] for j in range(i+1,N)),F(0)) for i in range(N)]
assert gamma(70)*max(visited_sumabs)<F(1,10)
assert gamma(2)*((1+gamma(70))*max(visited_sumabs)+2)<F(1,100)
assert max(visited)<2**40
assert N*N*int(np.max(np.abs(Q)))*max(visited)**2<2**127

# Exact dyadic denominator grids exclude subnormal intermediates.
Uden=max(x.denominator.bit_length()-1 for row in uf for x in row)
Dden=max(x.denominator.bit_length()-1 for x in df)
assert 2*Uden+Dden<1000
assert radius.denominator.bit_length()-1<1000

run=json.loads((P/'coset_full_4_logged.jsonl').read_text())
assert run['id']==4 and run['status']=='empty_float_candidate'
assert 0<run['nodes']<20_000_000_000 and run['witness']==[]

out={
 'status':'exact_global_error_bound_passed',
 'assumptions':['IEEE754 binary64 round-to-nearest operations; no unsafe fast-math','Actual logged coefficients and exact matrix are those supplied to the completed traversal','A full empty run must finish without node_cap'],
 'radius':enc(radius),'radius_margin':enc(margin),
 'max_coordinate_cap':max(M),'max_center_sum_abs':enc(max(s)),
 'max_stored_center_coefficient_error':enc(max(abs(uf[i][j]-U[i][j]) for i in range(N) for j in range(i+1,N))),
 'max_stored_diagonal_error':enc(max(abs(df[i]-d[i]) for i in range(N))),
 'max_center_error':enc(max(eta)),
 'component_error_sum':enc(component_sum),'partial_error':enc(partial_error),
 'max_branch_pruning_error':enc(branch_bound),'max_initial_base_pruning_error':enc(base_bound),
 'remaining_safety_margin':enc(margin-max(branch_bound,base_bound)),
 'all_visited_arithmetic_domain_max':max(visited),
 'all_visited_center_rounding_bound':enc(gamma(70)*max(visited_sumabs)),
 'dyadic_denominator_exponent_bound':2*Uden+Dden,
 'files_sha256':{name:hashlib.sha256((P/name).read_bytes()).hexdigest() for name in ['ne64_primal_best.npz','ne64_primal_best_Q.txt','ne64_certified_run_ldl.hex','coset_enum_logged.cpp','coset_full_4.txt']},
 'true_vector_coordinate_caps':M,'true_gso_component_caps':Z,
 'rounding_audit_is_not_traversal':True,
 'completed_logged_traversal_checked':True,
 'completed_logged_traversal_nodes':run['nodes'],
}
(P/'ne64_rounding_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['status','radius_margin','max_center_error','max_branch_pruning_error','max_initial_base_pruning_error','remaining_safety_margin','all_visited_arithmetic_domain_max']},indent=2))
