#!/usr/bin/env python3
"""Exact radial supporting planes for the unchanged accepted origin energy.

Every center is rationalized in integer-ray coordinates.  The proof is a
sum-of-squares completion followed by exact Bernstein restriction on each
physical radial piece.  This is a local floor certificate, not a global
density certificate.  It does not use a numerical SDP solution as evidence.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import json
import sys

if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT/'current/candidates/round9_six_graphs_degree18'
SCHEMA = 'RAY_CENTER_RADIAL_ORIGIN_FLOOR_V1'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rational(value):
    assert type(value) in (int, str), 'Exact inputs must be rational strings or integers.'
    return Q(value)


def split(values, t):
    assert Q(0) <= t <= Q(1)
    rows = [list(values)]
    while len(rows[-1]) > 1:
        rows.append([(1-t)*a+t*b for a, b in zip(rows[-1], rows[-1][1:])])
    return [row[0] for row in rows], [row[-1] for row in reversed(rows)]


def restrict(values, a, b):
    assert Q(0) <= a < b <= Q(1)
    _, right = split(values, a)
    left, _ = split(right, (b-a)/(1-a))
    return left


def nonnegative(values, bounds, node_limit):
    stack = [(list(values), bounds, 0)]
    nodes = leaves = maxdepth = 0
    while stack:
        controls, (a, b), depth = stack.pop()
        nodes += 1
        maxdepth = max(maxdepth, depth)
        assert nodes <= node_limit and depth <= 80, ('UNRESOLVED_RADIAL_POLYNOMIAL', a, b, depth)
        if min(controls) >= 0:
            leaves += 1
            continue
        assert controls[0] >= 0 and controls[-1] >= 0, ('NEGATIVE_EXACT_RADIAL_ENDPOINT', a, b)
        left, right = split(controls, Q(1, 2))
        assert left[-1] >= 0, ('NEGATIVE_EXACT_RADIAL_MIDPOINT', (a+b)/2, left[-1])
        stack.append((left, (a, (a+b)/2), depth+1))
        stack.append((right, ((a+b)/2, b), depth+1))
    return dict(status='EXACT_RADIAL_BERNSTEIN_NONNEGATIVE', nodes=nodes,
                leaves=leaves, maximum_depth=maxdepth)


def verify(record, source=SOURCE, node_limit=100000):
    assert __debug__, 'Run without -O.'
    source = Path(source).resolve()
    assert record['schema'] == SCHEMA
    kernel_path, origin_path = source/'kernel.json', source/'origin.json'
    assert record['kernel_sha256'] == sha(kernel_path)
    assert record['origin_sha256'] == sha(origin_path)
    kernel = json.loads(kernel_path.read_text())
    origin = json.loads(origin_path.read_text())
    assert kernel['dimension'] == 4 and kernel['local_degree'] == 6
    breaks = list(map(rational, kernel['breaks']))
    intervals = [list(map(rational, pair)) for pair in kernel['radial_basis_intervals']]
    assert len(breaks) == 8 and len(intervals) == 7 and breaks[0] == 1
    assert all(a < b for a, b in zip(breaks, breaks[1:]))
    assert intervals == [[Q(1), Q(3, 2)]]*6+[[Q(3, 2), Q(7, 4)]]
    assert breaks[-2:] == [Q(3, 2), Q(7, 4)]
    denominator = rational(kernel['common_denominator'])
    assert denominator > 0 and denominator.denominator == 1
    denominator = denominator.numerator
    rays = kernel['blocks'][0]['rays']
    centers = list(map(rational, record['ray_centers']))
    assert len(centers) == len(rays) and rays
    u = list(map(rational, origin['u']))
    v = list(map(rational, origin['v']))
    assert len(u) == len(v) == 49
    linear = [Q(0)]*49
    image_u = [Q(0)]*49
    constant = kappa = Q(0)
    for ray, center in zip(rays, centers):
        weight = ray['weight']
        vector = ray['vector']
        assert type(weight) is int and weight > 0
        assert len(vector) == 50 and vector[0] == 0
        assert all(type(x) is int for x in vector)
        vector = vector[1:]
        factor = Q(weight, denominator)
        vu = sum(x*y for x, y in zip(vector, u))
        constant += factor*(2*center*vu-center*center)
        kappa += factor*vu*vu
        for j in range(49):
            linear[j] += 2*factor*center*vector[j]
            image_u[j] += factor*vu*vector[j]
    assert image_u == v
    assert kappa == rational(origin['kappa']) == sum(a*b for a, b in zip(u, v))
    supplied = record.get('near_shell_lower_bounds')
    if supplied is not None:
        supplied = list(map(rational, supplied))
        assert len(supplied) == 5
    piece_records = []
    lower = []
    for piece, ((broad_a, broad_b), (a, b)) in enumerate(zip(intervals, zip(breaks, breaks[1:]))):
        assert broad_a <= a < b <= broad_b
        alpha, beta = (a-broad_a)/(broad_b-broad_a), (b-broad_a)/(broad_b-broad_a)
        controls = restrict(linear[7*piece:7*(piece+1)], alpha, beta)
        bound = (supplied[piece] if supplied is not None else min(controls)) if piece < 5 else Q(0)
        proof = nonnegative([x-bound for x in controls], (a, b), node_limit)
        if piece < 5:
            lower.append(bound)
        piece_records.append(dict(piece=piece, physical_interval=list(map(str, (a, b))),
                                  broad_basis_interval=list(map(str, (broad_a, broad_b))),
                                  exact_restricted_Bernstein_coefficients=list(map(str, controls)),
                                  certified_lower_bound=str(bound), proof=proof))
    count_linear = [lower[j]-lower[j+1] for j in range(4)]+[lower[4]]
    own = Path(__file__)
    binding = {str(path.relative_to(ROOT)): sha(path) for path in [kernel_path, origin_path, own]}
    result = dict(status='EXACT_RAY_CENTER_RADIAL_ORIGIN_FLOOR_VERIFIED',
                  global_density_bound_proved=False, schema=SCHEMA,
                  ray_centers=list(map(str, centers)), ray_count=len(rays),
                  physical_radial_linear_coefficients=list(map(str, linear)),
                  constant=str(constant), constant_decimal=float(constant),
                  near_shell_lower_bounds=list(map(str, lower)),
                  affine_count_coefficients=list(map(str, count_linear)),
                  exact_affine_floor_convention='constant + sum_i affine_count_coefficients[i]*N_i; N_i are the five cumulative near counts.',
                  residual_square_identity='E0-C-p*z=sum_i(weight_i/D)*(vector_i*(z+u)-center_i)^2>=0',
                  radial_sum_inequality='sum_x p(|x|)>=sum_j lower_j*(N_j-N_(j-1)); the two uncounted outer pieces have lower bound0.',
                  radial_endpoint_convention='Each physical piece is proved on its entire closed interval, including both one-sided values at every break.',
                  origin_range_equations=49, origin_square_identity_verified=True,
                  all_seven_closed_radial_pieces_verified=True, pieces=piece_records,
                  all_nonorigin_kernel_and_origin_couplings_unchanged=True,
                  source_sha256=binding)
    if 'target_counts' in record:
        counts = record['target_counts']
        assert len(counts) == 5 and all(type(x) is int and x >= 0 for x in counts)
        assert counts == sorted(counts)
        value = constant+sum(coefficient*N for coefficient, N in zip(count_linear, counts))
        result.update(target_counts=counts, target_pure_origin_floor=str(value),
                      target_pure_origin_floor_decimal=float(value))
    for path, digest in binding.items():
        assert sha(ROOT/path) == digest
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--source', type=Path, default=SOURCE)
    parser.add_argument('--node-limit', type=int, default=100000)
    args = parser.parse_args()
    assert not args.output.exists()
    proposal = json.loads(args.input.read_text())
    source = args.source.resolve()
    source_hashes = proposal.get('source_sha256', {})
    if source_hashes:
        for path, digest in source_hashes.items():
            assert sha(ROOT/path) == digest
        assert source_hashes[str((source/'kernel.json').relative_to(ROOT))] == sha(source/'kernel.json')
        assert source_hashes[str((source/'origin.json').relative_to(ROOT))] == sha(source/'origin.json')
    records = proposal.get('records', [proposal])
    exact_inputs = []
    reports = []
    for record in records:
        exact = dict(schema=SCHEMA, kernel_sha256=sha(source/'kernel.json'),
                     origin_sha256=sha(source/'origin.json'), ray_centers=record['ray_centers'])
        for optional in ('near_shell_lower_bounds', 'target_counts'):
            if optional in record:
                exact[optional] = record[optional]
        reports.append(verify(exact, source, args.node_limit))
        exact_inputs.append(exact)
    result = dict(status='EXACT_RAY_CENTER_RADIAL_ORIGIN_FLOOR_BUNDLE_VERIFIED',
                  global_density_bound_proved=False, exact_inputs=exact_inputs, records=reports,
                  proposal_file=str(args.input.resolve().relative_to(ROOT)), proposal_sha256=sha(args.input),
                  verifier_sha256=sha(Path(__file__)))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(status=result['status'], output=str(args.output), output_sha256=sha(args.output),
                         target_floors=[r.get('target_pure_origin_floor_decimal') for r in reports]), indent=2))


if __name__ == '__main__':
    main()
