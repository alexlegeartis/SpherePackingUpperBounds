#!/usr/bin/env python3
"""Exact maximum of origin lower bounds on the complete count domain.

The caller must freshly verify each radial plane, the projection equations,
count geometry and spherical surplus theorem before using this integer
enumeration in a global composition.  Standalone output is conditional.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import importlib.util
import json
import math
import sys

if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SOURCE = ROOT/'current/candidates/round9_six_graphs_degree18'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def enumerate_profile(projection, affine_records, spherical_profile, stratified_records=()):
    pure = list(map(Q, projection['projection_quadratic_coefficients']))
    total = list(map(Q, projection['projection_plus_existing_count_coefficients']))
    counts = [a-b for a, b in zip(total, pure)]
    assert len(pure) == len(total) == 21
    affine = [[Q(record['constant'])]+list(map(Q, record['affine_count_coefficients']))
              for record in affine_records]
    assert all(len(row) == 6 for row in affine)
    assert all(record['status'] == 'EXACT_RAY_CENTER_RADIAL_ORIGIN_FLOOR_VERIFIED'
               for record in affine_records)
    assert [row['N'] for row in spherical_profile] == list(range(27))
    spherical = [Q(row['spherical_surplus']) for row in spherical_profile]
    assert min(spherical) >= 0
    stratified = [(int(record['fixed_count']), int(record['radius_index']), Q(record['surplus']))
                  for record in stratified_records]
    assert all(0 <= N <= 26 and 0 <= j < 5 and surplus >= 0 for N, j, surplus in stratified)
    denominator = math.lcm(*(x.denominator for row in [pure, total, counts]+affine+[spherical]
                             for x in row), *(surplus.denominator for N, j, surplus in stratified))
    to_int = lambda row: [int(x*denominator) for x in row]
    total, count = to_int(total), to_int(counts)
    affine = [to_int(row) for row in affine]
    spherical = to_int(spherical)
    stratified = [(N, j, int(surplus*denominator)) for N, j, surplus in stratified]
    upper = [(i, j) for i in range(5) for j in range(i, 5)]
    profiles = {}
    base_profiles = {}
    first_profiles = {}
    checked = 0
    improved_vectors = 0

    def evaluate(coefficients, ns):
        return coefficients[0]+2*sum(coefficients[i+1]*ns[i] for i in range(5))+sum(
            coefficient*ns[i]*ns[j]*(1 if i == j else 2)
            for coefficient, (i, j) in zip(coefficients[6:], upper))

    def update(table, key, value, ns, winner, selected_surplus):
        if key not in table or (value, ns) < (table[key][0], table[key][1]):
            table[key] = (value, ns, winner, selected_surplus)

    for N in range(27):
        for n3 in range(min(N, 25)+1):
            for n2 in range(min(n3, 24)+1):
                for n1 in range(n2+1):
                    if n1 == 24 and n3 != 24:
                        continue
                    tail = (n1, n2, n3, N)
                    affine_constants = [row[0]+sum(a*n for a, n in zip(row[2:], tail)) for row in affine]
                    for n0 in range(n1+1):
                        if n0 == 24 and N != 24:
                            continue
                        ns = (n0,)+tail
                        old_total = evaluate(total, ns)
                        count_value = evaluate(count, ns)
                        supporting_values = [constant+row[1]*n0 for constant, row in zip(affine_constants, affine)]
                        if supporting_values:
                            winner = max(range(len(supporting_values)), key=supporting_values.__getitem__)
                            affine_total = count_value+supporting_values[winner]
                        else:
                            winner, affine_total = -1, old_total
                        if affine_total > old_total:
                            local_total = affine_total
                            improved_vectors += 1
                        else:
                            local_total = old_total
                            winner = -1
                        surplus = spherical[N]
                        for fixed_N, j, candidate in stratified:
                            if N == fixed_N and ns[j] == N:
                                surplus = max(surplus, candidate)
                        value = local_total+surplus
                        checked += 1
                        update(profiles, N, value, ns, winner, surplus)
                        update(base_profiles, N, local_total+spherical[N], ns, winner, spherical[N])
                        first = next(j+1 for j, n in enumerate(ns) if n == N)
                        update(first_profiles, (N, first), value, ns, winner, surplus)
    assert checked == projection['count_vectors_checked'] == 165503

    def rows(table):
        result = []
        for key, (value, ns, winner, surplus) in sorted(table.items()):
            row = dict(N=ns[-1], exact_floor=str(Q(value, denominator)),
                       floor_decimal=float(Q(value, denominator)), minimizing_counts=list(ns),
                       origin_lower_bound_winner='count_projection' if winner == -1 else f'radial_plane_{winner}',
                       spherical_surplus=str(Q(surplus, denominator)))
            if isinstance(key, tuple):
                row['first_containing_radius_index'] = key[1]
            result.append(row)
        return result

    minimum = min(value for value, ns, winner, surplus in profiles.values())
    before = min(value for value, ns, winner, surplus in base_profiles.values())
    return dict(status='EXACT_RADIAL_PROJECTION_MAXIMUM_COUNT_ENUMERATION_COMPLETE',
                global_density_bound_proved=False, count_vectors_checked=checked,
                affine_origin_floors=len(affine_records), improved_origin_floor_vectors=improved_vectors,
                origin_floor_rule='max(old exact E0 projection, all exact radial affine E0 floors)',
                existing_count_energies_added_exactly_once=True,
                spherical_rule='Maximum of applicable lower bounds for the same spherical energy.',
                conditional_count_profile=rows(profiles),
                first_containing_radius_profile=rows(first_profiles),
                profile_without_new_radius_stratification=rows(base_profiles),
                exact_origin_credit=str(Q(minimum, denominator)),
                exact_origin_credit_decimal=float(Q(minimum, denominator)),
                credit_without_new_radius_stratification=str(Q(before, denominator)),
                external_proof_dependencies_required=['Exact origin projection identities',
                    'Exact radial supporting-plane signs and integer ray completions',
                    'Geometric count bounds and implications',
                    'Source-matched spherical common-g surplus proofs with their exact radius/count hypotheses'])


if __name__ == '__main__':
    projection_path = HERE/'ACCEPTED_JOINT_COUNT_PROJECTION.json'
    radial_paths = [HERE/'RADIAL_FEASIBILITY_FLOORS_VERIFIED.json',
                    HERE/'RADIAL_FEASIBILITY_SUB32_FLOORS_VERIFIED.json']
    current_path = HERE/'accepted_uniform_joint_spherical_replay/GLOBAL_UNIFORM_JOINT_ORIGIN_CREDIT_VERIFIED.json'
    narrow_path = ROOT/'target625/geometry/accepted_stratified_N24_R1_face_exact_eps2/VERIFIED.json'
    current = json.loads(current_path.read_text())
    projection = json.loads(projection_path.read_text())
    local = load('fresh_envelope_projection', HERE/'verify_joint_count_origin_credit.py').verify(projection, SOURCE)
    assert local['count_vectors_checked'] == 165503
    radial_module = load('fresh_envelope_radial', HERE/'verify_ray_center_radial_floor.py')
    records = []
    for path in radial_paths:
        bundle = json.loads(path.read_text())
        for exact, prior in zip(bundle['exact_inputs'], bundle['records']):
            fresh = radial_module.verify(exact, SOURCE)
            assert fresh == prior
            records.append(fresh)
    narrow = json.loads(narrow_path.read_text())
    stratified = [dict(fixed_count=narrow['fixed_count'], radius_index=narrow['condition_radius_index'],
                       surplus=narrow['weighted_surplus_at_N'])]
    result = enumerate_profile(projection, records, current['conditional_count_profile'], stratified)
    raw = Q(current['exact_source_objective'])
    for key, credit in [('conditional_strict_threshold', Q(result['exact_origin_credit'])),
                        ('conditional_threshold_without_stratification', Q(result['credit_without_new_radius_stratification']))]:
        objective = raw-credit/8
        scaled = objective*10**12
        threshold = Q(scaled.numerator//scaled.denominator+1, 10**12)
        result[key] = str(threshold)
        result[key+'_decimal'] = float(threshold)
    paths = [Path(__file__), projection_path, current_path, narrow_path,
             HERE/'verify_joint_count_origin_credit.py', HERE/'verify_ray_center_radial_floor.py']+radial_paths
    result['source_sha256'] = {str(path.relative_to(ROOT)): sha(path) for path in paths}
    result['new_stratified_triangle_proof_freshly_replayed_here'] = False
    output = HERE/'RADIAL_ORIGIN_ENVELOPE_CONDITIONAL.json'
    assert not output.exists()
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(status=result['status'], credit=result['exact_origin_credit_decimal'],
                         conditional_threshold=result['conditional_strict_threshold_decimal'],
                         without_stratification=result['conditional_threshold_without_stratification_decimal'],
                         output_sha256=sha(output),
                         key_profiles=[row for row in result['conditional_count_profile'] if row['N'] in (24,25,26)]), indent=2))
