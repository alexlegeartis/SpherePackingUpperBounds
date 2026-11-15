"""Exact graph E with projected S2 degree-nine proof at beta=7/8.

The S3 annular count is replayed from its complete rational polynomial.
The S2 count uses an exactly factored quartic, with strictly positive
Legendre coefficients. No numerical neighbor count is assumed.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
import sys

from verify_exterior_graphs import parse, harmonic_power, prove_nonnegative, max_angle


def multiply(a, b):
    out = [Q(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def legendre_power(coefficients):
    # Normalized Legendre P_l, P_l(1)=1, from their exact recurrence.
    ps = [[Q(1)], [Q(0), Q(1)]]
    for degree in range(1, len(coefficients)-1):
        out = [Q(0)]*(degree+2)
        for i, value in enumerate(ps[-1]):
            out[i+1] += Q(2*degree+1, degree+1)*value
        for i, value in enumerate(ps[-2]):
            out[i] -= Q(degree, degree+1)*value
        ps.append(out)
    polynomial = [Q(0)]*len(coefficients)
    for coefficient, p in zip(coefficients, ps):
        for i, value in enumerate(p):
            polynomial[i] += coefficient*value
    return polynomial


def verify(data):
    assert data['schema'] == 'S3_PROJECTED_S2_ANTIPODAL_EXTERIOR_GRAPH_V1'
    assert data['global_density_bound_proved'] is False
    assert data['ambient_dimension'] == 4
    radial = tuple(map(parse, data['radial_interval']))
    angular = tuple(map(parse, data['target_angular_interval']))
    assert radial == (Q(33, 28), Q(7, 4))
    assert angular == (Q(-1), Q(-7, 8))
    L, R = radial
    rho = parse(data['spherical_cutoff'])
    assert rho == Q(41, 49)
    corners = [(r*r+s*s-1)/(2*r*s) for r in radial for s in radial]
    assert max(corners) == max_angle(radial, radial) == rho < 1

    # Fresh S3 count: a0*n^2 <= sum_ij p(ui.uj) <= n*p(1).
    coefficients = list(map(parse, data['count_harmonic_coefficients']))
    assert len(coefficients) >= 2
    assert coefficients[0] > 0 and min(coefficients[1:]) >= 0
    polynomial = harmonic_power(coefficients)
    sign = prove_nonnegative([-x for x in polynomial], Q(-1), rho)
    count_upper = sum(coefficients)/coefficients[0]
    count = count_upper.numerator // count_upper.denominator
    assert count == data['claimed_count_bound'] == 151

    beta = -angular[1]
    assert rho < beta < 1
    # For h_i,h_j >= beta, projected normalized dot is at most
    # F=(rho-h_i*h_j)/sqrt((1-h_i^2)(1-h_j^2)). Its partial derivative
    # has numerator rho*h_i-h_j <= rho-beta < 0. Hence F is maximized
    # at (beta,beta). If h_i=1, another neighbor is impossible because
    # it would have dot >= beta > rho with that axial neighbor.
    kappa = (rho-beta*beta)/(1-beta*beta)
    assert kappa == parse(data['projected_spherical_cutoff']) == Q(223, 735)
    assert 0 < kappa < 1
    legendre = list(map(parse, data['projected_legendre_coefficients']))
    assert len(legendre) == 5 and min(legendre) > 0
    shift = parse(data['projected_double_root_shift'])
    assert shift == Q(2, 3)
    factored = multiply(multiply([Q(1), Q(1)], [-kappa, Q(1)]),
                        multiply([shift, Q(1)], [shift, Q(1)]))
    assert legendre_power(legendre) == factored
    # On [-1,kappa], (t+1)>=0, (t-kappa)<=0, and (t+2/3)^2>=0.
    # Thus p<=0 throughout the full interval, with no root approximation.
    projected_upper = sum(legendre)/legendre[0]
    assert projected_upper == Q(2560, 279) < 10
    degree = projected_upper.numerator // projected_upper.denominator
    assert degree == data['claimed_degree_bound'] == 9
    unrounded_capacity = count*degree
    capacity = 2*(unrounded_capacity//2)
    assert unrounded_capacity == 1359
    assert capacity == data['claimed_even_ordered_edge_bound'] == 1358

    # Negative selected cosines make squared distance increasing in radii.
    distance_min = 2*L*L*(1+beta)
    assert distance_min == Q(16335, 3136) > R*R
    return dict(
        name=data['name'],
        status='EXACT_PROJECTED_S2_ANTIPODAL_EXTERIOR_GRAPH_CERTIFICATE_PASS',
        global_density_bound_proved=False,
        inherited_neighbor_count_used=False,
        radial_interval=data['radial_interval'],
        target_angular_interval=data['target_angular_interval'],
        minimum_selected_distance_squared=str(distance_min),
        minimum_exterior_distance_gap_squared=str(distance_min-R*R),
        spherical_cutoff=str(rho), radial_corner_cosines=list(map(str, corners)),
        count_proof=dict(
            method='S3 addition formula and exact whole-interval polynomial sign',
            f0=str(coefficients[0]), f1=str(sum(coefficients)),
            exact_count_upper=str(count_upper), count_bound=count,
            harmonic_coefficients=list(map(str, coefficients)),
            power_coefficients=list(map(str, polynomial)),
            whole_interval_sign_proof=sign,
        ),
        degree_proof=dict(
            cap_threshold=str(beta), projection_dimension=3,
            normalized_projected_pair_cosine_upper=str(kappa),
            monotonicity_margin_beta_minus_rho=str(beta-rho),
            axial_neighbor_case='If one neighbor projects to zero, the degree is at most1 because beta>rho.',
            projected_polynomial_factorization='(t+1)*(t-223/735)*(t+2/3)^2',
            projected_power_coefficients=list(map(str, factored)),
            projected_legendre_coefficients=list(map(str, legendre)),
            projected_f0=str(legendre[0]), projected_f1=str(sum(legendre)),
            projected_exact_count_upper=str(projected_upper),
            projected_gap_to_next_integer=str(degree+1-projected_upper),
            degree_bound=degree,
            whole_interval_sign_proof='Exact factorization; every factor sign proved throughout[-1,223/735].',
        ),
        unrounded_ordered_edge_bound=unrounded_capacity,
        even_ordered_edge_bound=capacity,
        parity_argument='Every undirected selected edge contributes exactly two ordered edges.',
        full_triangle_selected_anchor_count=0,
        exterior_selected_anchor_count_upper=1,
        kernel_diagonal='0', linear_term='0',
    )


if __name__ == '__main__':
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).with_name('antipodal_cap9_exterior_graph_certificate.json')
    report = verify(json.loads(path.read_text()))
    report['certificate_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    report['verifier_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    print(json.dumps(report, indent=2))
