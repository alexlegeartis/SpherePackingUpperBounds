"""Standalone rational proof of the near antipodal graph capacity 120.

The sign of the count polynomial follows from its explicit factorization.
All Gegenbauer identities, cutoff inequalities, integer bounds and support
comparisons are recomputed exactly.  This is a local theorem, not a global
packing certificate.  Only Python's standard library is used.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import json


def parse(value):
    assert isinstance(value, (str,int)) and not isinstance(value,bool)
    return Q(value)


def add(a,b):
    answer = [Q(0)]*max(len(a),len(b))
    for i,c in enumerate(a): answer[i] += c
    for i,c in enumerate(b): answer[i] += c
    return answer


def multiply(a,b):
    answer = [Q(0)]*(len(a)+len(b)-1)
    for i,c in enumerate(a):
        for j,d in enumerate(b): answer[i+j] += c*d
    return answer


def s3_zonals(degree):
    # U_0=1, U_1=2t, U_{k+1}=2t U_k-U_{k-1}; Z_k=U_k/(k+1).
    raw = [[Q(1)]]
    if degree: raw.append([Q(0),Q(2)])
    while len(raw) <= degree:
        raw.append(add([Q(0)]+[2*c for c in raw[-1]],[-c for c in raw[-2]]))
    result = [[c/Q(k+1) for c in row] for k,row in enumerate(raw)]
    assert all(sum(row)==1 for row in result)
    return result


def verify(data):
    assert data['schema'] == 'EXACT_FACTORED_S3_NEAR_ANTIPODAL_EXTERIOR_GRAPH_V1'
    assert data['dimension'] == 4
    lo,hi = map(parse,data['radial_interval'])
    tlo,thi = map(parse,data['target_angular_interval'])
    cutoff = parse(data['kernel_support_radius'])
    assert (lo,hi)==(Q(1),Q(103,100))
    assert (tlo,thi)==(Q(-1),Q(-3,4)) and cutoff==Q(7,4)
    beta = -thi
    # d/dr ((r^2+s^2-1)/(2rs)) has numerator r^2-s^2+1.
    # The lower bound below is positive on the whole radial square.
    derivative_numerator_lower = lo*lo-hi*hi+1
    assert derivative_numerator_lower > 0
    cosine = lambda r,s:(r*r+s*s-1)/(2*r*s)
    corner_cosines = [cosine(r,s) for r in (lo,hi) for s in (lo,hi)]
    rho = max(corner_cosines)
    assert rho==cosine(hi,hi)==parse(data['maximum_pair_cosine'])==Q(5609,10609)
    roots = list(map(parse,data['count_polynomial_roots']))
    assert roots==[rho,Q(-1),Q(-1,2),Q(-1,2),Q(0),Q(0)]
    polynomial = [Q(1)]
    for root in roots: polynomial = multiply(polynomial,[-root,Q(1)])
    coefficients = list(map(parse,data['normalized_s3_coefficients']))
    assert len(coefficients)==7 and all(c>0 for c in coefficients)
    expanded = [Q(0)]*len(polynomial)
    for coefficient,zonal in zip(coefficients,s3_zonals(6)):
        expanded = add(expanded,[coefficient*c for c in zonal])
    assert expanded==polynomial, 'The factored polynomial must match all zonal coefficients.'
    value_at_one = sum(polynomial)
    assert value_at_one==sum(coefficients)>0
    exact_count_upper = value_at_one/coefficients[0]
    assert exact_count_upper==parse(data['claimed_exact_count_upper'])==Q(288000,9391)
    count_bound = exact_count_upper.numerator//exact_count_upper.denominator
    assert count_bound==data['claimed_count_bound']==30
    assert 30<=exact_count_upper<31
    # On [-1,rho], (t-rho)<=0, (t+1)>=0, and both remaining
    # factors are squares. The S3 addition theorem therefore gives
    # p0*N^2 <= sum_ij p(u_i.u_j) <= N*p(1), including every N.
    projected_upper = rho-beta*beta
    assert projected_upper==parse(data['claimed_projected_pair_inner_product_upper'])==Q(-5737,169744)<0
    # For a fixed graph vertex u and selected neighbors v_i, put
    # h_i=-u.v_i>=beta and w_i=v_i+h_i*u in u-perp=R3.
    # Then w_i.w_j<=rho-h_i*h_j<=projected_upper<0.
    # Five vectors in R3 have an affine dependence with coefficients
    # summing to zero. Separate its positive and negative parts to get
    # equal positive combinations w. Their cross product is negative,
    # contradicting w.w>=0. Hence the degree is at most four.
    degree_bound = 4
    assert degree_bound==data['claimed_degree_bound']
    ordered = count_bound*degree_bound
    assert ordered==data['claimed_even_ordered_edge_bound']==120 and ordered%2==0
    minimum_distance_squared = 2*lo*lo*(1+beta)
    assert minimum_distance_squared==parse(data['claimed_minimum_selected_distance_squared'])==Q(7,2)
    assert minimum_distance_squared>cutoff*cutoff
    return dict(
        status='EXACT_NEAR_ANTIPODAL_EXTERIOR_GRAPH_CERTIFICATE_PASS',name=data['name'],
        radial_interval=list(map(str,(lo,hi))),target_angular_interval=list(map(str,(tlo,thi))),
        kernel_support_radius=str(cutoff),maximum_pair_cosine=str(rho),
        radial_geometry_proof=dict(corner_cosines=list(map(str,corner_cosines)),
            derivative_numerator_lower=str(derivative_numerator_lower),
            method='STRICT_MONOTONICITY_IN_EACH_RADIUS_ON_THE_ENTIRE_CLOSED_SQUARE'),
        count_proof=dict(count_bound=count_bound,exact_count_upper=str(exact_count_upper),
            polynomial_power_coefficients=list(map(str,polynomial)),
            normalized_s3_coefficients=list(map(str,coefficients)),
            constant_coefficient=str(coefficients[0]),value_at_one=str(value_at_one),
            coefficientwise_zonal_identity=True,
            whole_interval_sign='(t-rho)*(t+1)*(t+1/2)^2*t^2 <= 0 on the entire closed interval [-1,rho]',
            all_cardinalities=True,inherited_neighbor_count_used=False),
        degree_proof=dict(degree_bound=degree_bound,projection_dimension=3,
            projected_pair_inner_product_upper=str(projected_upper),
            method='AT_MOST_DIMENSION_PLUS_ONE_VECTORS_WITH_STRICTLY_NEGATIVE_MUTUAL_INNER_PRODUCTS'),
        exact_ordered_edge_upper=str(ordered),even_ordered_edge_bound=ordered,
        minimum_selected_distance_squared=str(minimum_distance_squared),
        distance_margin_over_kernel_support_squared=str(minimum_distance_squared-cutoff*cutoff),
        closed_radial_and_angular_endpoints=True,self_diagonal='0',linear_origin_coupling='0',
        energy='delta*(120-ordered_selected_edge_count)>=0 for every finite packing neighborhood',
        full_branch_decrement='0 because every selected edge is longer than the kernel support',
        global_density_bound_proved=False)


if __name__=='__main__':
    if not __debug__: raise RuntimeError('Run without -O.')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate',type=Path)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    certificate_bytes=args.certificate.read_bytes()
    report=verify(json.loads(certificate_bytes))
    report['certificate_sha256']=hashlib.sha256(certificate_bytes).hexdigest()
    report['verifier_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.output:
        assert not args.output.exists()
        args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
