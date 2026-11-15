"""Exact unrestricted single-color spherical graph bound in dimension four.

For P(t)=a0+sum_{k>=1} ak Zk(t), ak>=0, the spherical addition formula
gives a0*N*N <= sum_{i,j} P(u_i.u_j).  If P<=0 on all allowed distinct
inner products and P<=-1 on selected pairs, the ordered edge count obeys
e <= N*P(1)-a0*N*N.  The exact concave quadratic maximum over integers is
computed here.  No inherited bound on N and no cap projection is used.
"""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import hashlib
import json
import sys


def parse(value):
    assert isinstance(value,str)
    rational=Q(value)
    assert str(rational)==value
    return rational


def add(a,b):
    result=[Q(0)]*max(len(a),len(b))
    for i,value in enumerate(a):result[i]+=value
    for i,value in enumerate(b):result[i]+=value
    return result


def zonals(degree):
    assert type(degree) is int and degree>=0
    u=[[Q(1)],[Q(0),Q(2)]]
    for k in range(2,degree+1):
        u.append(add([Q(0)]+[2*x for x in u[-1]],[-x for x in u[-2]]))
    return [[x/Q(k+1) for x in polynomial] for k,polynomial in enumerate(u[:degree+1])]


def harmonic_power(coefficients):
    polynomial=[Q(0)]
    for coefficient,zonal in zip(coefficients,zonals(len(coefficients)-1)):
        polynomial=add(polynomial,[coefficient*x for x in zonal])
    return polynomial


def bernstein(polynomial,left,right):
    degree=len(polynomial)-1
    transformed=[Q(0)]*(degree+1)
    for j,value in enumerate(polynomial):
        for k in range(j+1):
            transformed[k]+=value*comb(j,k)*left**(j-k)*(right-left)**k
    return [sum((transformed[k]*Q(comb(i,k),comb(degree,k)) for k in range(i+1)),Q(0))
            for i in range(degree+1)]


def halves(coefficients):
    rows=[coefficients]
    while len(rows[-1])>1:
        rows.append([(left+right)/2 for left,right in zip(rows[-1],rows[-1][1:])])
    return [row[0] for row in rows],[row[-1] for row in reversed(rows)]


def prove_nonnegative(polynomial,left,right):
    assert left<=right
    stack=[(bernstein(polynomial,left,right),0)]
    nodes=leaves=maximum_depth=0
    while stack:
        coefficients,depth=stack.pop()
        nodes+=1;maximum_depth=max(maximum_depth,depth)
        assert nodes<=200000 and depth<80, 'Incomplete exact interval sign proof.'
        if min(coefficients)>=0:
            leaves+=1
            continue
        assert coefficients[0]>=0 and coefficients[-1]>=0, 'Negative exact endpoint.'
        assert max(coefficients)>=0, 'Negative entire interval.'
        left_child,right_child=halves(coefficients)
        stack.extend([(right_child,depth+1),(left_child,depth+1)])
    return dict(method='EXACT_RATIONAL_BERNSTEIN_SUBDIVISION',nodes=nodes,leaves=leaves,maximum_depth=maximum_depth)


def verify(data):
    assert data['schema']=='S3_SINGLE_COLOR_EXTERIOR_GRAPH_V1'
    assert data['global_density_bound_proved'] is False
    assert data['ambient_dimension']==4
    radial=tuple(map(parse,data['radial_interval']))
    angular=tuple(map(parse,data['target_angular_interval']))
    kernel_radius=parse(data['kernel_support_radius'])
    assert kernel_radius==Q(7,4)
    assert Q(1)<=radial[0]<radial[1]<=kernel_radius
    assert Q(-1)<=angular[0]<angular[1]<Q(0)
    assert parse(data['minimum_pair_distance_squared'])==1
    # For fixed s, the derivative numerator r^2-s^2+1 is increasing in r.
    # Any interior critical point is a minimum.  Applying this twice shows
    # the maximum over the radial square is attained at one of four corners.
    corner_cosines=[(r*r+s*s-1)/(2*r*s) for r in radial for s in radial]
    cutoff=max(corner_cosines)
    assert cutoff==parse(data['spherical_cutoff'])<1
    assert angular[1]<=cutoff
    coefficients=list(map(parse,data['harmonic_coefficients']))
    assert len(coefficients)==data['harmonic_degree']+1
    assert coefficients[0]>0 and min(coefficients[1:])>=0
    polynomial=harmonic_power(coefficients)
    assert polynomial==list(map(parse,data['power_coefficients']))
    # Zk(1)=1 follows directly from the stated recurrence.  This is also
    # checked against the actual constructed power polynomial.
    a0=coefficients[0];at_one=sum(coefficients,Q(0))
    assert sum(polynomial,Q(0))==at_one and at_one>0
    signs=[]
    for left,right,target in [(Q(-1),angular[0],0),
                              (angular[0],angular[1],1),
                              (angular[1],cutoff,0)]:
        nonnegative=[-x for x in polynomial]
        nonnegative[0]-=target
        signs.append(dict(left=str(left),right=str(right),target=target,
            proof=prove_nonnegative(nonnegative,left,right)))
    center=at_one/(2*a0)
    integer_center=center.numerator//center.denominator
    quadratic=lambda n: n*at_one-a0*n*n
    candidates=sorted({0,max(0,integer_center),max(0,integer_center+1)})
    maximum=max(quadratic(n) for n in candidates)
    capacity=2*(maximum.numerator//(2*maximum.denominator))
    assert capacity==data['claimed_even_ordered_edge_bound'] and capacity>=0
    # Negative selected cosine makes squared distance increasing in both
    # positive radii.  Thus the minimum occurs at the lower/lower corner.
    minimum_distance=2*radial[0]**2*(1-angular[1])
    assert minimum_distance==parse(data['minimum_selected_distance_squared'])>kernel_radius**2
    return dict(status='EXACT_SINGLE_COLOR_EXTERIOR_GRAPH_CERTIFICATE_PASS',
        name=data['name'],global_density_bound_proved=False,
        inherited_neighbor_count_used=False,ambient_dimension=4,
        radial_interval=data['radial_interval'],target_angular_interval=data['target_angular_interval'],
        kernel_support_radius=str(kernel_radius),spherical_cutoff=str(cutoff),
        radial_corner_cosines=list(map(str,corner_cosines)),
        harmonic_degree=data['harmonic_degree'],harmonic_coefficients=list(map(str,coefficients)),
        power_coefficients=list(map(str,polynomial)),
        whole_interval_sign_proofs=signs,
        edge_count_proof=dict(method='S3 addition formula and exact unrestricted concave quadratic maximum over nonnegative integer cardinalities',
            constant_coefficient=str(a0),diagonal_at_one=str(at_one),
            quadratic_center=str(center),real_quadratic_maximum=str(at_one*at_one/(4*a0)),
            integer_candidates=[dict(cardinality=n,upper_bound=str(quadratic(n))) for n in candidates],
            exact_integer_cardinality_edge_upper=str(maximum),
            gap_to_next_even=str(capacity+2-maximum),
            ordered_edge_evenness='Every unordered selected pair contributes its two opposite orders.'),
        even_ordered_edge_bound=capacity,
        minimum_selected_distance_squared=str(minimum_distance),
        minimum_exterior_distance_gap_squared=str(minimum_distance-kernel_radius**2),
        closed_radial_and_angular_endpoints=True,
        full_triangle_selected_anchor_count=0,exterior_selected_anchor_count_upper=1,
        kernel_diagonal='0',linear_term='0')


if __name__=='__main__':
    path=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('moderate_exterior_graph_certificate.json')
    report=verify(json.loads(path.read_text()))
    report['certificate_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
    report['verifier_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    print(json.dumps(report,indent=2))
