#!/usr/bin/env python3
"""Fresh global composition using the maximum of common-g surplus proofs.

The frozen fixed-N composer first replays the accepted 111+1 global proof,
the count projection, its geometric dependencies and both fixed-N proofs.
This wrapper then freshly checks the bounded-N26 proof for the same energy
and takes its pointwise maximum with the fixed-N lower bounds.  The new
total origin credit replaces the old credit once; pair/cyclic data are fixed.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import sys
import time

if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
GEO = ROOT / 'target625/geometry'
SOURCE = ROOT / 'current/candidates/round9_six_graphs_degree18'
FIXED = HERE / 'compose_joint_count_global_credit.py'
CHECKER = GEO / 'certify_accepted_combined_boundedN_credit.py'
CERT = GEO / 'accepted_combined_boundedN26_exact_eps7_v2/certificate.json'
PRIOR = CERT.with_name('VERIFIED.json')
ALGEBRA = ROOT / 'optimization/continuation_20261018/base/new_work/spherical_sdp/certify.py'
AUDITOR = HERE / 'audit_bounded_cardinality_repair.py'
PINS = {
    FIXED: '4db7677b8d0ef9d242a324c36af2ff4f06d3e518894391e1ce4b0bb4889f2960',
    CHECKER: '9296ccbfd4468de4f2738b78b78d6625eced56e5ef5a2d4f46498c9e81e75cfb',
    CERT: '50f5632f6262683069bc0fcf7507feda8deb2a4ef6c3b63f901350dc415a78c9',
    PRIOR: 'ebfb65eac32a725ee189eb9eb9b3a477f4cbf5d2b9b6a69efee8e2baaa15088b',
    ALGEBRA: 'aeb873edff74112cd0dd0293d1c77be91aec7ec405aa1013b9b92088db46fca5',
    AUDITOR: 'aebadb20c2c74d36200afed9623457f7dd5c90ab83dba4dd5eb4d37dc8f63b50',
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write(path, value):
    temporary = path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(value, indent=2)+'\n')
    temporary.replace(path)


def stable(value):
    if isinstance(value, dict):
        return {key: stable(item) for key, item in value.items() if key != 'seconds'}
    if isinstance(value, list):
        return [stable(item) for item in value]
    return value


def compose(output_dir):
    assert __debug__, 'Run without -O.'
    started = time.time()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=False)
    binding = {str(path.relative_to(ROOT)): digest for path, digest in PINS.items()}
    for path, digest in binding.items():
        assert sha(ROOT/path) == digest, path
    own = Path(__file__)
    binding[str(own.relative_to(ROOT))] = sha(own)
    # This is an actual fresh accepting composition, with all original gates.
    fixed_module = load('fixed_joint_count_composer_fresh', FIXED)
    inherited, inherited_path = fixed_module.compose(output_dir/'fixed_replay')
    assert inherited['status'] == 'EXACT_GLOBAL_FIVE_COUNT_AND_FIXED_CARDINALITY_ORIGIN_CREDIT_VERIFIED'
    assert inherited['global_density_bound_proved'] is True
    assert inherited['original142_payloads_unchanged'] is True
    assert inherited['complete_domain_coverage']['total_branches'] == 112
    assert inherited['count_geometry_freshly_replayed'] is True
    assert inherited['fixedN_proofs_freshly_replayed'] == 2
    assert inherited['all_pair_and_cyclic_expressions_unchanged'] is True
    assert inherited['all_graphs_retained_unchanged'] is True
    for path, digest in inherited['proof_file_sha256'].items():
        assert sha(ROOT/path) == digest, path
        assert path not in binding or binding[path] == digest, path
        binding[path] = digest
    binding[str(inherited_path.relative_to(ROOT))] = sha(inherited_path)
    structural = load('independent_bounded_count_audit_fresh', AUDITOR).verify()
    assert structural['status'] == 'INDEPENDENT_EXACT_BOUNDED_N26_COMMON_G_REPAIR_AUDIT_PASS'
    assert structural['maximum_count'] == 26
    structural_path = output_dir/'FRESH_BOUNDED_N26_STRUCTURAL_AUDIT.json'
    write(structural_path, structural)
    binding[str(structural_path.relative_to(ROOT))] = sha(structural_path)
    for path, digest in structural['source_sha256'].items():
        assert path not in binding or binding[path] == digest, path
        binding[path] = digest
    bounded = load('accepted_bounded_count_fresh', CHECKER)
    assert Path(bounded.algebra.__file__).resolve() == ALGEBRA.resolve()
    assert sha(Path(bounded.algebra.__file__)) == PINS[ALGEBRA]
    assert sha(Path(bounded.common_record.__code__.co_filename)) == fixed_module.FIXED_PINS[
        GEO/'certify_accepted_combined_fixedN_credit.py']
    fresh = bounded.verify(CERT, node_limit=1200000)
    assert fresh['status'] == 'EXACT_BOUNDED_CARDINALITY_SURPLUS_OF_ACCEPTED_SPHERICAL_SUM_VERIFIED'
    assert fresh['maximum_code_count'] == 26 and fresh['global_bound_proved'] is False
    assert fresh['pair']['status'] == fresh['triangle']['status'] == 'CERTIFIED'
    assert fresh['common_g_unchanged'] is True
    assert fresh['certificate_sha256'] == PINS[CERT]
    assert fresh['verifier_sha256'] == PINS[CHECKER]
    assert fresh['algebra_sha256'] == PINS[ALGEBRA]
    assert stable(fresh) == stable(json.loads(PRIOR.read_text()))
    origin = json.loads((SOURCE/'origin.json').read_text())
    rays = origin['appended_spherical_local_energy_rays']
    assert len(rays) == 2
    degree = len(fresh['common_scalar_coefficients'])
    assert degree == 13
    scale = sum(Q(ray['weight']) for ray in rays)
    g = [sum(Q(ray['weight'])*Q(ray['scalar_coefficients'][j]) for ray in rays)
         for j in range(degree)]
    B = sum(Q(ray['weight'])*Q(ray['B']) for ray in rays)
    combination = fresh['accepted_combination']
    assert combination['accepted_origin_sha256'] == sha(SOURCE/'origin.json')
    assert list(map(Q, combination['weighted_accepted_g'])) == g
    assert Q(combination['weighted_accepted_B']) == B
    assert Q(combination['positive_normalization']) == scale > 0
    assert list(map(Q, fresh['common_scalar_coefficients'])) == [value/scale for value in g]
    per_neighbor = B-scale*Q(fresh['normalized_B'])
    assert per_neighbor == Q(fresh['weighted_surplus_per_neighbor']) > 0
    assert per_neighbor == Q(structural['weighted_surplus_per_neighbor'])
    fresh_path = output_dir/'FRESH_BOUNDED_N26_SURPLUS.json'
    write(fresh_path, fresh)
    binding[str(fresh_path.relative_to(ROOT))] = sha(fresh_path)
    old_profile = inherited['conditional_count_profile']
    assert [row['N'] for row in old_profile] == list(range(27))
    assert min(Q(row['total_floor']) for row in old_profile) == Q(inherited['new_origin_credit'])
    profile = []
    for row in old_profile:
        N = row['N']
        base = Q(row['base_count_projection_floor'])
        fixed = Q(row['spherical_surplus'])
        assert Q(row['total_floor']) == base+fixed and fixed >= 0
        assert N in (24, 25) or fixed == 0
        uniform = N*per_neighbor
        # All candidates bound the SAME unchanged spherical energy.
        surplus = max(Q(0), uniform, fixed)
        winner = 'fixed_N' if fixed >= uniform and fixed > 0 else 'bounded_N26'
        profile.append(dict(N=N, base_count_projection_floor=str(base),
                            fixed_cardinality_surplus=str(fixed),
                            bounded_cardinality_surplus=str(uniform),
                            spherical_surplus=str(surplus), chosen_lower_bound=winner,
                            total_floor=str(base+surplus),
                            minimizing_nested_counts=row['minimizing_nested_counts']))
    credit = min(Q(row['total_floor']) for row in profile)
    raw = Q(inherited['exact_source_objective'])
    fixed_credit = Q(inherited['new_origin_credit'])
    original_credit = Q(inherited['old_origin_credit'])
    constant = Q(inherited['graph_origin_constant_retained'])
    assert raw == Q(origin['objective'])+constant/8
    assert raw == Q(inherited['objective'])+fixed_credit/8
    assert credit > fixed_credit > original_credit > 0
    objective = raw-credit/8
    scaled = objective*10**12
    threshold = Q(scaled.numerator//scaled.denominator+1, 10**12)
    assert objective < threshold < Q(inherited['strict_upper_bound'])
    local = dict(status='EXACT_FIVE_COUNT_AND_MAX_COMMON_G_LOCAL_CREDIT_VERIFIED',
                 global_density_bound_proved=False, exact_origin_credit=str(credit),
                 exact_density_credit=str(credit/8),
                 count_vectors_checked=165503, complete_count_profile=profile,
                 five_count_projection_and_geometry_certificate=inherited['local_certificate'],
                 fresh_fixed_N_proofs=2, fresh_bounded_N26_proof=True,
                 bounded_N26_surplus_per_neighbor=str(per_neighbor),
                 spherical_lower_bound='max(0,N*d,certified fixed-N24/N25 surplus)',
                 same_energy_estimates_summed=False,
                 aggregate_local_energy='E0 + the three unchanged count energies + the same two unchanged spherical energies - new credit >= 0.',
                 residual_square_identity='E0-(N+m)^T D_inverse (N+m)=(z+u-B D_inverse (N+m))^T M (z+u-B D_inverse (N+m)) >= 0',
                 old_credit_accounting='This new total replaces every previous origin credit exactly once.',
                 all_pair_and_cyclic_expressions_unchanged=True)
    local_path = output_dir/'LOCAL_MAX_COMMON_G_ORIGIN_CREDIT_VERIFIED.json'
    write(local_path, local)
    binding[str(local_path.relative_to(ROOT))] = sha(local_path)
    for path, digest in binding.items():
        assert sha(ROOT/path) == digest, path
    record = dict(status='EXACT_GLOBAL_FIVE_COUNT_AND_MAX_COMMON_G_ORIGIN_CREDIT_VERIFIED',
                  global_density_bound_proved=True, source=str(SOURCE.relative_to(ROOT)),
                  inherited_fresh_global_proof=str(inherited_path.relative_to(ROOT)),
                  inherited_fresh_global_proof_sha256=sha(inherited_path),
                  original_accepted_global_proof_sha256=inherited['original_accepted_global_proof_sha256'],
                  original142_payloads_unchanged=True,
                  complete_domain_coverage=inherited['complete_domain_coverage'],
                  all_pair_and_cyclic_expressions_unchanged=True,
                  all_graphs_retained_unchanged=True, nonorigin_kernel_unchanged=True,
                  total_physical_origin_coupling_unchanged=True,
                  graph_origin_constant_retained=str(constant),
                  graph_density_cost_retained=inherited['graph_density_cost_retained'],
                  old_origin_credit=str(original_credit), replaced_fixed_only_credit=str(fixed_credit),
                  new_origin_credit=str(credit), exact_density_credit=str(credit/8),
                  additional_density_improvement=str((credit-original_credit)/8),
                  improvement_over_fixed_only_density=str((credit-fixed_credit)/8),
                  exact_source_objective=str(raw), objective=str(objective),
                  objective_decimal=float(objective), strict_upper_bound=str(threshold),
                  threshold_decimal=float(threshold), conditional_count_profile=profile,
                  local_certificate=str(local_path.relative_to(ROOT)),
                  local_certificate_sha256=sha(local_path),
                  count_geometry_freshly_replayed=True, fixedN_proofs_freshly_replayed=2,
                  boundedN26_proof_freshly_replayed=True,
                  same_spherical_energy_lower_bounds_combined_by='pointwise maximum',
                  origin_credit_accounting='Replace the old total credit, without adding either competing spherical lower bounds or previous origin credits.',
                  theorem='The freshly verified count projection and common-g surpluses give an all-cardinality nonnegative local decomposition with a smaller independent origin constant. All complete inherited pair/cyclic and Fourier conditions are unchanged, so the strict density bound applies to arbitrary four-dimensional sphere packings.',
                  new_global_interval_run_required=False, new_global_interval_run_claimed=False,
                  proof_file_sha256=binding, composer_source=str(own.relative_to(ROOT)),
                  composer_sha256=sha(own), elapsed_seconds=time.time()-started)
    output = output_dir/'GLOBAL_UNIFORM_JOINT_ORIGIN_CREDIT_VERIFIED.json'
    write(output, record)
    return record, output


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    record, path = compose(args.output_dir)
    print(json.dumps(dict(status=record['status'], strict_upper_bound=record['strict_upper_bound'],
                         objective_decimal=record['objective_decimal'],
                         output=str(path), output_sha256=sha(path)), indent=2))
