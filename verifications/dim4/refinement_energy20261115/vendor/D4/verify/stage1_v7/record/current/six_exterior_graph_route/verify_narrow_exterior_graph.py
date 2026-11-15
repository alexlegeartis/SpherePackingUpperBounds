"""Exact graph E certificate: a fresh S3 count plus a degree-four cap lemma.

No retained neighbor count is assumed.  Rational Bernstein subdivision proves
the complete spherical polynomial sign condition.  The only geometric lemma
used afterward is the following elementary fact: at most d+1 vectors in R^d
can have pairwise strictly negative inner products.  If there were d+2, an
affine dependence would split into two nonempty positive combinations with
the same sum w.  Their cross inner product would be strictly negative, whereas
it equals ||w||^2.  Here d=3, after projection perpendicular to a graph vertex.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
import sys

from verify_exterior_graphs import (
    parse, harmonic_power, prove_nonnegative, max_angle,
)


def verify(data):
    assert data['schema'] == 'S3_ANTIPODAL_EXTERIOR_GRAPH_V2'
    assert data['global_density_bound_proved'] is False
    assert data['ambient_dimension'] == 4
    radial = tuple(map(parse, data['radial_interval']))
    angular = tuple(map(parse, data['target_angular_interval']))
    assert radial == (Q(33, 28), Q(7, 4))
    assert angular == (Q(-1), Q(-183, 200))
    L, R = radial
    cutoff = parse(data['spherical_cutoff'])
    assert cutoff == Q(41, 49)
    corners = [(r*r+s*s-1)/(2*r*s) for r in radial for s in radial]
    assert max(corners) == max_angle(radial, radial) == cutoff < 1

    # Addition formula: a0*n^2 <= sum_{i,j} p(u_i.u_j) <= n*p(1).
    # This proves the count for the actual broad annulus, from scratch.
    coefficients = list(map(parse, data['count_harmonic_coefficients']))
    assert len(coefficients) >= 2
    assert coefficients[0] > 0 and min(coefficients[1:]) >= 0
    polynomial = harmonic_power(coefficients)
    sign = prove_nonnegative([-x for x in polynomial], Q(-1), cutoff)
    f0, f1 = coefficients[0], sum(coefficients)
    count_upper = f1/f0
    count = count_upper.numerator // count_upper.denominator
    assert count == data['claimed_count_bound'] == 151

    # For p=-u and every selected neighbor v, v.p >= a.  Distinct selected
    # neighbors have projected inner product <= cutoff-a^2 < 0 in p^perp.
    cap_threshold = -angular[1]
    projected_upper = cutoff-cap_threshold*cap_threshold
    assert projected_upper == -Q(961, 1960000) < 0
    projection_dimension = data['ambient_dimension']-1
    degree = projection_dimension+1
    assert degree == data['claimed_degree_bound'] == 4
    capacity = count*degree
    assert capacity == data['claimed_even_ordered_edge_bound'] == 604
    assert capacity % 2 == 0

    # Every selected edge is outside the supported full-triangle branch.
    distance_min = 2*L*L*(1-angular[1])
    assert distance_min == Q(417087, 78400) > R*R
    return dict(
        name=data['name'],
        status='EXACT_ANTIPODAL_EXTERIOR_GRAPH_CERTIFICATE_PASS',
        global_density_bound_proved=False,
        inherited_neighbor_count_used=False,
        radial_interval=data['radial_interval'],
        target_angular_interval=data['target_angular_interval'],
        minimum_selected_distance_squared=str(distance_min),
        minimum_exterior_distance_gap_squared=str(distance_min-R*R),
        spherical_cutoff=str(cutoff),
        radial_corner_cosines=list(map(str, corners)),
        count_proof=dict(
            method='S3 addition formula and exact whole-interval polynomial sign',
            f0=str(f0), f1=str(f1), exact_count_upper=str(count_upper),
            count_bound=count, gap_to_next_integer=str(count+1-count_upper),
            harmonic_coefficients=list(map(str, coefficients)),
            power_coefficients=list(map(str, polynomial)),
            whole_interval_sign_proof=sign,
        ),
        degree_proof=dict(
            cap_threshold=str(cap_threshold),
            projection_dimension=projection_dimension,
            projected_pair_inner_product_upper=str(projected_upper),
            degree_bound=degree,
            lemma='Pairwise strictly negative inner products in R^d imply cardinality at most d+1, by affine dependence.',
        ),
        even_ordered_edge_bound=capacity,
        full_triangle_selected_anchor_count=0,
        exterior_selected_anchor_count_upper=1,
        kernel_diagonal='0', linear_term='0',
    )


if __name__ == '__main__':
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).with_name('antipodal_wide_exterior_graph_certificate.json')
    report = verify(json.loads(path.read_text()))
    report['certificate_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    report['verifier_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    print(json.dumps(report, indent=2))
