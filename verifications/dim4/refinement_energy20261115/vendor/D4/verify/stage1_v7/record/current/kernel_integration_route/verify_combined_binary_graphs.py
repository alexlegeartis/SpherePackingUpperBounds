"""Exact verification using the same complete coordinator and a compiled
integer midpoint subdivision. The C++20 recurrence is entrywise identical
to the original NumPy recurrence; it changes no polynomial, bound, sign
decision, subdivision axis, or coverage requirement. The backend is built
fresh from the included source. Source and binary are both bound before
and after the run. Angular degree reduction is a separate data transform.
"""
"""Exact verification with independent antipodal allocation and a count gauge.

This is a separate coordinator cloned from the immutable audited sequential
verify_count_gauge.py.  The classical and two-level antipodal corrections use a separate exact helper.
The pair-origin identity includes its full linear cost before the complete
count budget is constructed. Earlier verifiers are unchanged. Scheduling
and every all-branch coverage gate are inherited from the audited parallel
coordinator.  The candidate's local exact_bernstein_fast.py is
hashed before and after the run along with all other proof inputs.

Before any fork, the parent freshly verifies the geometric/count dependencies,
the scalar, the exact origin and objective, and every complete base pair.
Each worker receives only (kind, regions), inherits the already validated
modules and budget, and verifies six independent polynomial/Bernstein
identities for its branch before calling the unchanged cyclic checker.
The exact crosschecks use a contracted transform. Every coefficient
used for signs is separately enclosed by directed integer arithmetic. A
positive rational binary scale is retained even when its exponent is negative.
Every expected branch must be present exactly once and CERTIFIED.  For seven
regions this means 84 full and 28 exterior branches and 672 identities.
The pair and cyclic stages are explicitly partial and never write VERIFIED.

All nonorigin counts, including whether the probe contains the distinguished
origin, are checked by the general ball-count contract.  Pointed corrections
are retained on every full and exterior triangle branch and on the pair self
diagonal.  The shifted PSD origin, the restored linear coupling, count
polynomial and exact objective are independently recomputed.

Affine latitude exclusion is replayed on its whole two-variable domain and
depends on a freshly replayed spherical count proof.  Negative-pole tents
use exact capacity bounds.  Their linear origin rows and all mixed terms
are recomputed here; certificate prose is never used as a formula.

This variant also verifies a full-pair count gauge: every base pair row is
freshly proved, the entire residual is retained, and its zero-cost local
energy cancels every supported pair expression exactly.

This separate variant selects the widest current exact coordinate interval
for subdivision. On exterior branches the third coordinate is normalized z.
All enclosure, acceptance, source, and complete-coverage gates are unchanged.
Both the immutable parent coordinator and selected helper are explicitly
pinned below and retained in the beginning/end source contract.

No numerical solver output, failed branch or incomplete stage is acceptance.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import concurrent.futures
import hashlib
import importlib.util
import itertools
import json
import math
import multiprocessing
import os
import platform
import sys
import time

INTEGRATION = Path(__file__).resolve().parent
BASE = next(p for p in [INTEGRATION.parent/'base']+[q/'historical_snapshot' for q in INTEGRATION.parents]
            if (p/'new_work/spherical_sdp/certify.py').is_file())
if (INTEGRATION.parent/'deps').is_dir():
    sys.path.insert(0, str(INTEGRATION.parent/'deps'))
sys.path.insert(0, str(BASE/'new_work/spherical_sdp'))
sys.path.insert(0, str(BASE/'new_work/pointed'))
import numpy as np
from count_quadratic import NestedCountCone
import signed_spherical_exact as signed_helper

HERE = BASE/'new_work/spherical_sdp'
_WORKER_STATE = None


def expected_branch_jobs(number_of_regions):
    """The same exhaustive unordered branch set as the sequential checker."""
    assert isinstance(number_of_regions, int) and number_of_regions >= 1
    jobs = tuple((kind, tuple(regions))
                 for kind, size in [('full', 3), ('outside', 2)]
                 for regions in itertools.combinations_with_replacement(
                     range(number_of_regions), size))
    full = math.comb(number_of_regions+2, 3)
    outside = math.comb(number_of_regions+1, 2)
    assert len(jobs) == len(set(jobs)) == full+outside
    if number_of_regions == 7:
        assert (full, outside, len(jobs), 6*len(jobs)) == (84, 28, 112, 672)
    return jobs


def validate_branch_packet(job, packet, expected, seen, angular_degree):
    """Fail closed on missing, duplicate, mislabelled or incomplete results."""
    key = packet['kind'], tuple(packet['regions'])
    assert key == job and key in expected and key not in seen, \
        'A branch result was duplicated, unexpected, or assigned to another job.'
    identity = packet['identities']
    assert identity['status'] == 'VERIFIED'
    assert identity['polynomial_identity_checks'] == 6
    assert identity['angular_degree'] == angular_degree
    assert identity['verification_scope'] == 'EXACT_ARITHMETIC_CROSSCHECKS_ONLY_NOT_A_GLOBAL_CERTIFICATE'
    assert identity['basis_transform_method'] == 'EXACT_TRANSPOSE_AXIS_CONTRACTION'
    assert identity['interval_converter_implementation_checked'] is False
    cyclic = packet['cyclic']
    signed_identity = cyclic['signed_spherical_identity']
    assert signed_identity['status'] == 'EXACT_SIGNED_SPHERICAL_BRANCH_IDENTITIES_VERIFIED'
    assert signed_identity['kind'] == key[0] and tuple(signed_identity['regions']) == key[1]
    assert signed_identity['checks'] == 5
    assert cyclic['kind'] == key[0] and tuple(cyclic['regions']) == key[1]
    assert cyclic['status'] == 'CERTIFIED', cyclic
    enclosure = cyclic['bernstein_enclosure']
    assert enclosure['method'] == 'EXACT_OUTWARD_TWO_PRODUCT_BERNSTEIN_V3'
    assert enclosure['all_coefficients_enclosed'] is True
    assert enclosure['coefficient_precision'] == 192
    assert isinstance(enclosure['output_exponent'], int)
    assert -4096 <= enclosure['output_exponent'] <= 50
    assert F(enclosure['exact_positive_output_scale']) == F(2)**enclosure['output_exponent'] > 0
    seen.add(key)
    return key


def assert_complete_branch_coverage(jobs, packets, angular_degree):
    """Replay the coverage gate independently of asynchronous arrival order."""
    expected, seen = set(jobs), set()
    assert len(jobs) == len(expected) and len(packets) == len(jobs)
    identity_count = 0
    for packet in packets:
        job = packet['kind'], tuple(packet['regions'])
        validate_branch_packet(job, packet, expected, seen, angular_degree)
        identity_count += packet['identities']['polynomial_identity_checks']
    assert seen == expected and identity_count == 6*len(jobs)
    return identity_count


def _branch_worker(job):
    """Use inherited modules; no polynomial tensor or module is pickled."""
    state = _WORKER_STATE
    assert state is not None and state['parent_ready']
    assert job in state['expected_jobs']
    kind, regions = job
    core, budget = state['core'], state['budget']
    assert budget.core is core and budget.nonnegative_proven
    started = time.time()
    try:
        identities = state['identities'].verify(state['folder'], only=job)
        assert identities['status'] == 'VERIFIED'
        assert identities['polynomial_identity_checks'] == 6
        assert identities['angular_degree'] == state['angular_degree']

        cone, coefficients = state['cone'], state['coefficients']
        active_flags = [state['flags'](region) for region in regions]
        if kind == 'outside':
            active_flags.append([0]*cone.k)
        correction = dot(cone.triangle_row_from_flags(active_flags), coefficients)
        joint_cyclic = F(0)
        for polynomial, weight in state['joint_polynomials']:
            joint_flags = [polynomial.flags(core.BREAKS, region) for region in regions]
            if kind == 'outside':
                joint_flags.append([0]*polynomial.k)
            joint_cyclic += weight*dot(polynomial.triangle_row(joint_flags),
                                      polynomial.coefficients)
        correction += joint_cyclic
        mixed = conditional_correction(core, kind, regions, state['conditional'])
        eta = -correction-mixed
        record = state['lift_helper'].check(
            core, kind, regions, eta, state['tents'], state['basis_intervals'],
            budget, state['lift_probes'], caps=state['caps'],
            allocations=state['allocations'], spherical_rays=state['spherical_rays'],
            graph_weights=state['graph_weights'])
        # UNKNOWN can be returned before the old helper attaches branch tags.
        # Attach our job tags, but reject any conflicting tags from the helper.
        assert record.get('kind', kind) == kind
        assert tuple(record.get('regions', regions)) == regions
        record.update(kind=kind, regions=list(regions), eta=str(eta),
                      pointed_correction='EXACT_CAP_TENT_CONDITIONAL_AND_TWO_CAPACITY_ANTIPODAL_BOUNDS',
                      conditional_polynomial_cyclic_correction=str(joint_cyclic),
                      conditional_cyclic_correction=str(mixed),
                      nested_cyclic_correction=str(correction), count_flags=active_flags)
        return dict(kind=kind, regions=list(regions), identities=identities,
                    cyclic=record, worker_pid=os.getpid(), seconds=time.time()-started)
    finally:
        # Retain reuse within a branch, but do not accumulate its large exact
        # tensors across jobs. Clearing a cache changes no proof computation.
        core._base_build.cache_clear()
        if hasattr(core._exact_backend, 'clear_cache'):
            core._exact_backend.clear_cache()


def _atomic_json(path, value):
    temporary = path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(value, indent=2)+'\n')
    temporary.replace(path)


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def add(a, b):
    result = [F(0)]*max(len(a), len(b))
    for i, value in enumerate(a):
        result[i] += value
    for i, value in enumerate(b):
        result[i] += value
    return result


def dot(a, b):
    return sum((F(x)*y for x, y in zip(a, b)), F(0))


def conditional_flags(core, region, ray):
    """Inner and open-inner/closed-outer shell flags on one physical piece.

    Both one-sided values at every jump are covered by the branch loop.
    For the actual local energy use closed inner balls, hence shell=N_R-N_r.
    """
    lo, hi = core.BREAKS[region:region+2]
    r, R = ray['inner_radius'], ray['outer_radius']
    assert hi <= r or lo >= r
    assert hi <= R or lo >= R
    return int(hi <= r), int(lo >= r and hi <= R)


def conditional_correction(core, kind, regions, rays):
    """The exact nonpositive mixed kernel, summed on a triangle."""
    assert kind in ['full', 'outside']
    assert len(regions) == (3 if kind == 'full' else 2)
    pairs = [(0, 1), (0, 2), (1, 2)] if kind == 'full' else [(0, 1)]
    result = F(0)
    for ray in rays:
        flags = [conditional_flags(core, region, ray) for region in regions]
        result -= ray['weight']/2*sum(flags[i][0]*flags[j][1]+flags[i][1]*flags[j][0]
                                     for i, j in pairs)
    return result


def validate_conditional(records, core, affine):
    """Validate the geometry and collect, without assuming, spherical proofs.

    A full verification later supplies freshly checked spherical reports.
    Pair/cyclic probes only hash this dependency and cannot certify a bound.
    """
    result, dependencies, reports = [], [], []
    for record in records:
        assert record['proof_kind'] == 'affine_latitude'
        assert record.get('kind', 'one_outer') == 'one_outer'
        assert record.get('maximum', 24) == 24
        assert F(record.get('pair_coefficient', 23)) == 23
        assert F(record.get('mixed_kernel_coefficient', '-1/2')) == F(-1, 2)
        r, R, beta = map(F, [record['inner_radius'], record['outer_radius'], record['weight']])
        # A positive-width inner piece is required to cover its closed
        # boundary value at contact; r=1 would need a separate singleton row.
        assert 1 < r < R <= core.BREAKS[-1] and R*R < 2 and beta >= 0
        assert r in core.BREAKS and R in core.BREAKS
        proof = record['proof']
        assert F(proof['inner_radius']) == r and F(proof['outer_radius']) == R
        assert proof['maximum_inner_count'] == 24
        checked = affine.check(r, R, F(proof['spherical_cutoff']))
        assert checked['status'] == 'EXACT_AFFINE_LATITUDE_DOMAIN_CERTIFIED'
        for key in ['inner_cosine_bound', 'pole_cosine_bound', 'beta', 'alpha',
                    'transformed_pole_upper']:
            assert F(checked[key]) == F(proof[key]), (key, checked[key], proof[key])
        assert F(checked['transformed_pole_upper']) <= F(proof['spherical_cutoff'])
        dependency = {'certificate': proof['spherical_certificate'],
                      'certificate_sha256': proof['spherical_certificate_sha256']}
        if 'spherical_certificate_sha256' in record:
            assert record['spherical_certificate_sha256'] == dependency['certificate_sha256']
        path = HERE/dependency['certificate']
        assert path.parent.resolve() == HERE.resolve()
        assert hashlib.sha256(path.read_bytes()).hexdigest() == dependency['certificate_sha256']
        result.append(dict(inner_radius=r, outer_radius=R, weight=beta,
                           proof=proof, dependency=dependency))
        dependencies.append(dependency)
        reports.append(dict(geometry=checked, spherical_dependency=dependency,
                            energy='(23-N_inner)*(N_outer-N_inner)',
                            weight=str(beta)))
    return result, dependencies, reports


def verify(folder, stage='all', workers=1):
    global _WORKER_STATE
    if not __debug__:
        raise RuntimeError('Run without -O.')
    assert stage in ['all', 'pair', 'cyclic']
    assert isinstance(workers, int) and 1 <= workers <= 4
    assert _WORKER_STATE is None, 'This coordinator cannot run overlapping verifications.'
    folder = Path(folder).resolve()
    started = time.time()
    input_names = ['origin.json', 'scalar.json', 'kernel.json', 'graph_energies.json']
    code_names = ['certificate_core.py', 'verify_scalar.py', 'audit_identities.py',
                  'exact_bernstein_fast.py']
    input_hashes = {name: hashlib.sha256((folder/name).read_bytes()).hexdigest()
                    for name in input_names+code_names}
    split_backend = load_module('compiled_integer_midpoint', INTEGRATION.parent/'fast_split.py')
    split_contract = split_backend.source_contract()
    contract_paths = [INTEGRATION.parent/'fast_split.py', INTEGRATION.parent/'fast_split.cpp',
                      split_backend.LIBRARY, Path(__file__), INTEGRATION/'signed_spherical_exact.py',
                      INTEGRATION/'verify_combined_binary_fastsplit.py',
                      INTEGRATION.parent/'origin_route/interval_combined_binary_graphs.py',
                      INTEGRATION.parent/'graph_route/verify_exterior_graphs.py',
                      INTEGRATION.parent/'graph_route/exterior_graph_certificates.json',
                      INTEGRATION.parent/'graph_route/graph_repair_bounds.py',
                      INTEGRATION/'interval_combined.py',
                      INTEGRATION/'interval_combined_refined.py',
                      INTEGRATION/'interval_combined_binary.py',
                      INTEGRATION/'verify_combined_binary.py',
                      INTEGRATION.parent/'origin_route/interval_combined_binary_physical_width.py',
                      INTEGRATION/'verify_combined_refined.py',
                      INTEGRATION.parent/'geometry_route/fixed_binary_bernstein_v3.py',
                      INTEGRATION.parent/'geometry_route/fixed_binary_bernstein_two_product_v3.py',
                      INTEGRATION.parent/'geometry_route/binary_bernstein_signed_exponent_v3.py',
                      HERE.parent/'pointed/exact_bernstein_fast.py',
                      INTEGRATION.parent/'origin_route/contracted_bernstein_identity_v2.py',
                      INTEGRATION.parent/'origin_route/audit_identities_contracted_v2.py',
                      INTEGRATION/'verify_combined.py',
                      INTEGRATION.parent/'origin_route/cancellation_preserving_pair_budget.py',
                      HERE/'verify_count_gauge.py',
                      HERE/'verify_count_gauge_fast_parallel.py',
                      HERE/'verify_independent_antipodal.py',
                      HERE.parent/'pointed'/'independent_antipodal_allocation_exact.py',
                      HERE.parent/'pointed'/'interval_certificate_independent_antipodal.py',
                      HERE.parent/'pointed'/'two_capacity_antipodal_allocation_exact.py',
                      HERE.parent/'pointed'/'interval_certificate_two_capacity_antipodal.py',
                      HERE/'count_quadratic.py', HERE/'certify.py',
                      HERE.parent/'pointed'/'interval_certificate.py',
                      HERE.parent/'pointed'/'interval_certificate_general_closed.py',
                      HERE.parent/'pointed'/'interval_certificate_full_closed.py',
                      HERE.parent/'pointed'/'pair_interval_certificate_closed.py',
                      HERE.parent/'pointed'/'interval_certificate_tent_closed.py',
                      HERE.parent/'pointed'/'pair_slack_lift.py',
                      HERE.parent/'pointed'/'interval_certificate_lift.py',
                      HERE.parent/'pointed'/'full_pair_count_gauge.py',
                      HERE.parent/'pointed'/'interval_certificate_count_gauge.py',
                      HERE.parent/'conditional'/'affine_latitude_counts.py',
                      HERE.parent/'conditional'/'joint_count_polynomial.py']
    contract_hashes = {str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest()
                       for path in contract_paths}
    for pinned_path, expected_hash in (
            (INTEGRATION/'verify_combined_binary.py',
             '1cca3a3dc908e57853e8b03c79d1ebec448a241567196024018264c5ad81dc14'),
            (INTEGRATION.parent/'origin_route/interval_combined_binary_physical_width.py',
             '15de750f695531141321e87d5e198e1baf48aaa71bd933916551c3e418bbbd96')):
        assert contract_hashes[str(pinned_path.resolve())] == expected_hash, \
            'The selected scheduling helper or immutable parent source changed.'
    sys.path.insert(0, str(folder))
    core = load_module('certificate_core', folder/'certificate_core.py')
    scalar = load_module('verify_scalar', folder/'verify_scalar.py')
    # Identical directed integer midpoint recurrence, using the source-bound compiled backend.
    core.split = split_backend.split
    assert core.ROOT.resolve() == folder
    assert Path(core._exact_backend.__file__).resolve() == folder/'exact_bernstein_fast.py'
    assert core.transform is core._exact_backend.transform
    assert core.to_bernstein is core._exact_backend.to_bernstein
    import flint
    runtime_versions = dict(python=platform.python_version(), numpy=np.__version__,
                            python_flint=flint.__version__)
    runtime_versions['midpoint_subdivision'] = split_contract
    jobs = expected_branch_jobs(core.NS)
    data = json.loads((folder/'origin.json').read_text())
    S = json.loads((folder/'scalar.json').read_text())
    # The scalar proof must cover the entire complement of kernel support.
    # A larger scalar-tail cutoff would leave an unverified radial gap.
    assert 1 <= F(S['cutoff']) <= core.BREAKS[-1]
    kernel_data = json.loads((folder/'kernel.json').read_text())
    graph_data = json.loads((folder/'graph_energies.json').read_text())
    assert graph_data['schema'] == 'INDEPENDENT_EXTERIOR_GRAPH_ENERGIES_V1'
    graph_weights = tuple(map(F, graph_data['weights']))
    assert len(graph_weights) == 3 and min(graph_weights) >= 0 and graph_weights[2] == 0
    graph_dir = INTEGRATION.parent/'graph_route'
    for name, digest in graph_data['source_sha256'].items():
        assert name in ['verify_exterior_graphs.py','exterior_graph_certificates.json','graph_repair_bounds.py']
        assert hashlib.sha256((graph_dir/name).read_bytes()).hexdigest() == digest
    assert set(graph_data['source_sha256']) == {'verify_exterior_graphs.py','exterior_graph_certificates.json','graph_repair_bounds.py'}
    graph_bounds = load_module('graph_repair_bounds', graph_dir/'graph_repair_bounds.py')
    graph_verifier = load_module('exact_exterior_graph_verifier', graph_dir/'verify_exterior_graphs.py')
    graph_certificates = json.loads((graph_dir/'exterior_graph_certificates.json').read_text())['certificates']
    assert len(graph_certificates) == 2
    graph_reports = [graph_verifier.verify(d) for d in graph_certificates]
    assert [r['even_ordered_edge_bound'] for r in graph_reports] == [820,1102]
    assert graph_bounds.CAPACITIES == (820,1102,9630)
    assert graph_bounds.NEAR == (F(1),F(103,100))
    assert graph_bounds.MIDDLE == (F(103,100),F(3,2))
    assert graph_bounds.OUTER == (F(3,2),F(7,4))
    for report,far,angular in zip(graph_reports,[graph_bounds.MIDDLE,graph_bounds.OUTER],
                                  [(F(-3,4),F(-1,2)),(F(-19,20),F(-3,4))]):
        assert report['status'] == 'EXACT_LOCAL_GRAPH_CERTIFICATE_PASS'
        assert tuple(tuple(map(F,row)) for row in report['radial_closure_intervals']) == (graph_bounds.NEAR,far)
        assert tuple(map(F,report['target_angular_interval'])) == angular
        assert F(report['minimum_selected_distance_squared']) > core.BREAKS[-1]**2
    assert core.BREAKS[-1] == F(7,4)
    graph_origin_constant = sum((w*r['even_ordered_edge_bound'] for w,r in zip(graph_weights,graph_reports)),F(0))
    assert graph_origin_constant/8 == graph_bounds.density_cost(graph_weights)
    graph_energy_report = dict(status='FRESH_EXACT_INDEPENDENT_GRAPH_ENERGIES_VERIFIED',
        weights=list(map(str,graph_weights)), proofs=graph_reports,
        origin_constant=str(graph_origin_constant), self_diagonal='0', linear_origin_coupling='0',
        energy='sum_j delta_j*(C_j - ordered_selected_edge_count_j)',
        full_branch_decrement='0: every selected edge is longer than the kernel support',
        pair_and_full_pair_gauge_unchanged=True)
    print('EXACT_INDEPENDENT_GRAPH_ENERGIES_PASSED', graph_energy_report, flush=True)
    nested = data['nested_count_quadratic']
    cone = NestedCountCone(nested['radii'], nested['bounds'])
    coefficients = list(map(F, nested['coefficients']))
    count_polynomial = cone.verify(coefficients)
    assert all(radius in core.BREAKS for radius in cone.radii)
    assert all(F(1) <= radius < core.BREAKS[-1] for radius in cone.radii)
    proof_records = nested['count_proofs']
    assert len(proof_records) == cone.k
    for radius, maximum, proof in zip(cone.radii, cone.bounds, proof_records):
        assert radius <= F(proof['R']) and proof['M'] <= maximum
    cuts = data.get('count_cuts', [])
    assert all(F(cut['weight']) >= 0 for cut in cuts)
    # Additional unverified energy representations must not be silently
    # ignored by a verifier implementing only this particular contract.
    assert 'integer_quadratic' not in data
    assert 'local_energy' not in data
    assert data['spherical_local_energy'] == 'SIGNED_SPHERICAL_REGISTRY_V1'
    assert not data.get('count_gauge_products', [])
    assert not data.get('product_count_gauges', [])
    sys.path.insert(0, str(HERE.parent/'pointed'))
    pointed_path = HERE.parent/'pointed'/'interval_certificate_full_closed.py'
    pointed = load_module('pointed_general_full_exact', pointed_path)
    pair_pointed = load_module('pointed_general_pair_exact', HERE.parent/'pointed'/'pair_interval_certificate_closed.py')
    tent_helper = load_module('geometric_tent_exact', HERE.parent/'pointed'/'interval_certificate_tent_closed.py')
    budget_helper = load_module('full_pair_count_gauge_exact', HERE.parent/'pointed'/'full_pair_count_gauge.py')
    refined_budget_helper = load_module('cancellation_preserving_pair_budget_exact',
        INTEGRATION.parent/'origin_route/cancellation_preserving_pair_budget.py')
    allocation_helper = load_module('two_capacity_antipodal_allocation_exact',
        HERE.parent/'pointed'/'two_capacity_antipodal_allocation_exact.py')
    lift_helper = load_module('widest_coordinate_binary_combined_signed_spherical_cyclic_exact',
        INTEGRATION.parent/'origin_route/interval_combined_binary_graphs.py')
    assert lift_helper.cyclic_decrement_lower is graph_bounds.cyclic_decrement_lower
    assert lift_helper.allocation_bounds is allocation_helper.correction_bounds
    affine = load_module('geometric_affine_exact', HERE.parent/'conditional'/'affine_latitude_counts.py')
    joint_helper = load_module('geometric_joint_count_exact', HERE.parent/'conditional'/'joint_count_polynomial.py')
    assert data['pointed_caps'], 'This contract requires explicit pointed corrections.'
    basis_intervals = [tuple(map(F, row))
                       for row in kernel_data['radial_basis_intervals']]
    caps = pointed.validate_caps(data['pointed_caps'], core, basis_intervals)
    tents = tent_helper.validate_tents(data.get('pointed_tents', []), core, basis_intervals)
    allocation_records = data['independent_antipodal_allocations']
    assert isinstance(allocation_records, list) and allocation_records
    allocations = allocation_helper.validate_allocations(
        allocation_records, core, basis_intervals)
    assert len(allocations) == len(allocation_records)
    spherical_local = load_module('signed_local_spherical_certifier', HERE/'certify.py')
    spherical_rays, spherical_ray_reports, signed_dependencies = signed_helper.validate_registry(
        data.get('appended_spherical_local_energy_rays', []), core, folder, spherical_local, replay=True)
    assert lift_helper.build_signed is signed_helper.build_signed
    assert not data.get('tail_budget_caps', []), 'This variant spends each pair slack only once.'
    gauge_record = data['full_pair_count_gauge']
    lift_probes = []
    conditional, conditional_dependencies, conditional_reports = validate_conditional(
        data.get('conditional_radial_rays', []), core, affine)
    joint_polynomials, joint_reports, joint_count_proofs = [], [], []
    all_geometry, all_geometry_reports = list(conditional), list(conditional_reports)
    for record in data.get('conditional_count_polynomials', []):
        polynomial = joint_helper.ConditionalCountPolynomial(record)
        weight = F(record['weight'])
        assert weight >= 0
        assert all(1 < r < core.BREAKS[-1] and r in core.BREAKS for r in polynomial.radii)
        count_proofs = record['count_proofs']
        assert len(count_proofs) == polynomial.k
        for radius, maximum, proof in zip(polynomial.radii, polynomial.bounds, count_proofs):
            assert radius <= F(proof['R']) and proof['M'] <= maximum
        joint_count_proofs.extend(count_proofs)
        specs = [dict(inner_radius=str(polynomial.radii[row['inner_index']]),
                      outer_radius=str(polynomial.radii[row['outer_index']]),
                      maximum=row['maximum'], weight='0', proof_kind=row['proof_kind'],
                      proof=row['proof']) for row in record['implications']]
        implications, dependencies, geometry_reports = validate_conditional(specs, core, affine)
        for report in geometry_reports:
            report['role'] = 'INTEGER_DOMAIN_IMPLICATION_ONLY'
        all_geometry.extend(implications)
        all_geometry_reports.extend(geometry_reports)
        conditional_dependencies.extend(dependencies)
        joint_polynomials.append((polynomial, weight))
        joint_reports.append(dict(integer_polynomial=polynomial.verify(),
                                  geometric_implications=geometry_reports, weight=str(weight),
                                  origin_constant=str(weight*polynomial.coefficients[0])))
    gauge_radii = list(cone.radii)
    gauge_bounds = list(cone.bounds)
    for polynomial, _ in joint_polynomials:
        gauge_radii.extend(polynomial.radii)
        gauge_bounds.extend(polynomial.bounds)
    gauge = budget_helper.validate_gauge(gauge_record, core, gauge_radii, gauge_bounds)
    ordinary_dependencies = proof_records+cuts+joint_count_proofs
    all_dependencies = ordinary_dependencies+conditional_dependencies
    for proof in all_dependencies:
        path = HERE/proof['certificate']
        assert path.parent.resolve() == HERE.resolve()
        assert hashlib.sha256(path.read_bytes()).hexdigest() == proof['certificate_sha256']
    result = {'status': 'INCOMPLETE_EXACT_FULL_PAIR_WIDEST_COORDINATE_BINARY_COMBINED_SIGNED_SPHERICAL_COUNT_GAUGE_FAST_PARALLEL_VERIFICATION',
              'count_polynomial': count_polynomial, 'pairs': [], 'cyclic': [],
              'signed_spherical_local_energies': spherical_ray_reports,
              'independent_graph_local_energies': graph_energy_report,
              'pointed_cap_count': len(caps),
              'pointed_tent_count': len(tents),
              'independent_antipodal_allocation_count': len(allocations),
              'independent_antipodal_capacity_proofs': [a['capacity_proof'] for a in allocations],
              'independent_antipodal_origin_constant': '0',
              'independent_antipodal_self_diagonal': '0',
              'full_pair_count_gauge': {key: [str(v) for v in value] if isinstance(value, list) else str(value)
                                        for key, value in gauge.items()},
              'tent_capacity_lemmas': [tent['capacity_lemma'] for tent in tents],
              'conditional_rays': conditional_reports,
              'conditional_count_polynomials': joint_reports,
              'truncated_pointed_cap_count': sum('pivot_lower_cutoff' in cap for cap in caps),
              'pointed_contract': 'GENERAL_SIGMA_AND_ACTIVE_NONORIGIN_BALL_COUNTS',
              'pointed_helper_sha256': hashlib.sha256(pointed_path.read_bytes()).hexdigest(),
              'initial_input_sha256': input_hashes,
              'contract_sha256': contract_hashes,
              'runtime_versions': runtime_versions,
              'exact_bernstein_fast_sha256': input_hashes['exact_bernstein_fast.py'],
              'parallel': dict(workers=workers, start_method='fork',
                               number_of_regions=core.NS,
                               expected_full_branches=math.comb(core.NS+2, 3),
                               expected_outside_branches=math.comb(core.NS+1, 2),
                               expected_jobs=[dict(kind=kind, regions=list(regions))
                                              for kind, regions in jobs],
                               expected_identity_checks=6*len(jobs),
                               accepted_jobs=0, branch_packets=[])}
    progress_path = folder/('WIDEST_COORDINATE_BINARY_COMBINED_SIGNED_SPHERICAL_COUNT_GAUGE_FAST_PARALLEL_PROGRESS.json' if stage == 'all' else
                            'WIDEST_COORDINATE_BINARY_COMBINED_SIGNED_SPHERICAL_COUNT_GAUGE_FAST_PARALLEL_'+stage.upper()+'_PROBE.json')

    def checkpoint():
        result['seconds'] = time.time()-started
        _atomic_json(progress_path, result)

    def assert_sources_unchanged():
        assert all(signed_helper.sha(path) == digest for path, digest in signed_dependencies), \
            'A signed spherical source changed during verification.'
        assert all(hashlib.sha256((folder/name).read_bytes()).hexdigest() == value
                   for name, value in input_hashes.items()), \
            'Certificate inputs or the candidate exact backend changed during verification.'
        assert all(hashlib.sha256(Path(path).read_bytes()).hexdigest() == value
                   for path, value in contract_hashes.items()), \
            'Verifier sources changed during verification.'
        assert all(hashlib.sha256((HERE/proof['certificate']).read_bytes()).hexdigest()
                   == proof['certificate_sha256'] for proof in all_dependencies), \
            'A spherical count certificate changed during verification.'

    checkpoint()
    if stage in ['all', 'cyclic']:
        spherical = load_module('spherical_nested_exact', HERE/'certify.py')
        checked, proof_results = {}, []
        # Each copy is byte-identical; replay outputs cannot overwrite older
        # accepted reports.  The originals are hashed again at completion.
        replay_folder = folder/'spherical_replay_fast_parallel'
        replay_folder.mkdir(exist_ok=True)
        for proof in all_dependencies:
            name = proof['certificate']
            if name not in checked:
                replay = replay_folder/name
                replay.write_bytes((HERE/name).read_bytes())
                checked[name] = spherical.verify(replay)
                assert checked[name]['certificate_sha256'] == proof['certificate_sha256']
        for proof in ordinary_dependencies:
            name = proof['certificate']
            report = checked[name]
            radius = F(proof['R'])
            assert 1 <= radius and radius*radius-radius-1 < 0
            assert 1-F(1)/(2*radius*radius) <= F(report['cosine_cutoff'])
            assert report['count_bound'] <= proof['M']
            proof_results.append({'certificate': name,
                                  'certificate_sha256': proof['certificate_sha256'],
                                  'R': proof['R'], 'M': proof['M'], 'verification': report})
        result['count_proofs'] = proof_results
        for ray, report in zip(all_geometry, all_geometry_reports):
            spherical_report = checked[ray['dependency']['certificate']]
            assert spherical_report['status'] == 'EXACT_SPHERICAL_COUNT_24_VERIFIED'
            assert spherical_report['count_bound'] <= 24
            assert F(ray['proof']['spherical_cutoff']) <= F(spherical_report['cosine_cutoff'])
            assert spherical_report['certificate_sha256'] == ray['dependency']['certificate_sha256']
            report['spherical_verification'] = spherical_report
            report['status'] = 'EXACT_AFFINE_CONDITIONAL_ENERGY_VERIFIED'
        result['scalar'] = scalar.verify()
        result['parent_scalar_and_dependencies'] = 'FRESH_EXACT_REPLAY_PASSED'
        # Each branch worker invokes this unchanged audit with an exact
        # (kind, regions) selector, before checking its inequality.
        # This audit retains all six exact evaluations and harmonic/
        # radial endpoint checks. Coefficient enclosure is separately supplied
        # by the pinned outward converter, not inferred from these evaluations.
        identities = load_module('audit_identities_contracted_v2',
            INTEGRATION.parent/'origin_route/audit_identities_contracted_v2.py')
        checkpoint()

    matrices, denominator, angular_degree = core.gram()
    v, u = [list(map(F, data[key])) for key in ['v', 'u']]
    kappa = F(data['kappa'])
    n = len(v)
    assert len(u) == n == core.NC*core.NS
    assert all(sum((F(int(matrices[0][i+1, j+1]), denominator)*u[j]
                    for j in range(n)), F(0)) == v[i] for i in range(n))
    assert dot(v, u) == kappa >= 0
    result['origin_identity'] = 'EXACT_MU_EQUALS_V_AND_VU_EQUALS_KAPPA'
    pair_v = list(map(F, data['pair_origin_v']))
    assert len(pair_v) == n
    origin_shift = [sum((cap['M']*cap['weight']*cap['w'][i]/2
                         for cap in caps), F(0)) for i in range(n)]
    for i in range(n):
        origin_shift[i] += sum((tent['capacity']*tent['weight']*tent['w'][i]/2
                               for tent in tents), F(0))
        origin_shift[i] += sum((a['weight']*a['w'][i]/2 for a in allocations), F(0))
        region = i//core.NC
        origin_shift[i] += sum((ray['weight']*ray['B']*ray['flags'][region]/2
                                for ray in spherical_rays), F(0))
        origin_shift[i] += sum((F(23, 2)*ray['weight']*conditional_flags(core, region, ray)[1]
                               for ray in conditional), F(0))
    assert all(v[i]+origin_shift[i] == pair_v[i] for i in range(n)), \
        'The auxiliary pair row must equal the PSD row plus every exact geometric linear coefficient.'
    result['pointed_origin_cancellation'] = 'EXACT_CAP_TENT_CONDITIONAL_AND_TWO_CAPACITY_ANTIPODAL_LINEAR_ORIGIN_IDENTITIES'
    pair_polynomials = {}
    verified_pair_majorants = []
    budget = None
    # Both nonorigin arguments have the same annulus support.  The auxiliary
    # row uses untruncated linear weights.  The pair reduction retains each
    # self diagonal, and removes the whole M*w term below a validated pivot
    # cutoff.  Thus it represents the actual (possibly truncated) energy.

    def flags(region):
        lo, hi = core.BREAKS[region:region+2]
        out = [int(hi <= radius) for radius in cone.radii]
        assert all(hi <= radius or lo >= radius for radius in cone.radii)
        return out

    fs, alpha = list(map(F, S['F'])), F(S['alpha'])
    objective = (fs[0]+kappa+coefficients[0]+
                 sum((cut['M']*F(cut['weight']) for cut in cuts), F(0))+
                 sum((weight*polynomial.coefficients[0]
                      for polynomial, weight in joint_polynomials), F(0)))/8
    assert alpha == 2 and F(S['G'][0]) == 1
    assert objective == F(data['objective']) < F(data['threshold'])
    graph_objective = objective+graph_origin_constant/8
    assert graph_objective < F(graph_data['threshold'])
    result.update(base_objective=str(objective), base_strict_upper_bound=data['threshold'],
                  objective=str(graph_objective), objective_decimal=float(graph_objective),
                  strict_upper_bound=graph_data['threshold'])
    if stage in ['all', 'pair'] or gauge is not None:
        for region in range(core.NS):
            left, right = core.BREAKS[region:region+2]
            boundaries = sorted(set([left, right]+[F(cut['R']) for cut in cuts
                                                 if left < F(cut['R']) < right]))
            B = core.radial(region)
            j = core.NC*region
            coupling = [sum((pair_v[j+i]*int(B[i, k]) for i in range(core.NC)), F(0))
                        for k in range(core.DEG+1)]
            diagonal = [F(0)]*(2*core.DEG+1)
            bound = 2*max(abs(x) for x in pair_v[j:j+core.NC])
            for ell, matrix in enumerate(matrices):
                M = matrix[1:, 1:] if ell == 0 else matrix
                M = M[j:j+core.NC, j:j+core.NC]
                bound += max(F(int(M[i, i]), denominator) for i in range(core.NC))
                polynomial = B.T@M@B
                for i, k in itertools.product(range(core.NC), repeat=2):
                    diagonal[i+k] += F(int(polynomial[i, k]), denominator)
            correction = dot(cone.pair_row_from_flags(flags(region)), coefficients)
            joint_pair = sum((weight*dot(polynomial.pair_row(polynomial.flags(core.BREAKS, region)),
                                        polynomial.coefficients)
                              for polynomial, weight in joint_polynomials), F(0))
            correction += joint_pair
            # B*chi is already present through the exact origin identity.
            # Retain the signed self term -g(1)*chi separately, exactly once.
            spherical_diagonal = -sum((ray['weight']*sum(ray['g'])*ray['flags'][region]
                                      for ray in spherical_rays), F(0))
            correction += spherical_diagonal
            assert sum((ray['weight']*(ray['B']-sum(ray['g']))*ray['flags'][region]
                        for ray in spherical_rays), F(0)) == sum(
                            (ray['weight']*ray['kappa']*ray['flags'][region]
                             for ray in spherical_rays), F(0))
            uncut_D = add(diagonal, [2*x for x in coupling])
            uncut_D[0] += correction
            pair_polynomials[region] = uncut_D
            for lo, hi in zip(boundaries, boundaries[1:]):
                mid = (lo+hi)/2
                active = sum((F(cut['weight']) for cut in cuts if mid <= F(cut['R'])), F(0))
                D = list(uncut_D)
                D[0] -= active
                positive_bound = bound+abs(correction)+active
                Fpoly = [F(0)]*(2*len(fs)-1)
                for k, coefficient in enumerate(fs):
                    Fpoly[2*k] = coefficient*alpha**k
                order = 55
                exponential = [F(0)]*(2*order+1)
                for k in range(order+1):
                    exponential[2*k] = alpha**k/math.factorial(k)
                upper = alpha*hi**2
                assert 0 < upper < order+2
                remainder = [F(0)]*(2*(order+1)+1)
                remainder[-1] = positive_bound*alpha**(order+1)/(math.factorial(order+1)*(1-upper/(order+2)))
                polynomial = add(add(Fpoly, core.conv(D, exponential)), remainder)
                pd = math.lcm(*(x.denominator for x in polynomial))
                pc = np.array([int(x*pd) for x in polynomial], dtype=object)
                pair_majorant_coefficients = tuple(map(int, pc))
                record = pair_pointed.check(core, pc, pd, (lo, hi), region, caps,
                                             basis_intervals, alpha)
                record.update(region=region, left=str(lo), right=str(hi),
                              nested_pair_correction=str(correction),
                              conditional_polynomial_pair_correction=str(joint_pair),
                              signed_spherical_diagonal=str(spherical_diagonal),
                              count_flags=flags(region), active_linear_weight=str(active),
                              exponential_order=order, positive_part_bound=str(positive_bound))
                result['pairs'].append(record)
                print('GEOMETRIC_POINTED_PAIR', record, flush=True)
                checkpoint()
                assert record['status'] == 'CERTIFIED', record
                # These are the actual coefficients passed immediately above
                # to the fresh pair verifier. No saved polynomial/status loads.
                assert tuple(map(int, pc)) == pair_majorant_coefficients
                verified_pair_majorants.append(dict(
                    region=region, left=lo, right=hi,
                    coefficients=pair_majorant_coefficients, denominator=pd))

    # Even a cyclic-only probe must freshly prove every complete unlifted
    # pair branch. A status claimed in an input file cannot enable the lift.
    budget = budget_helper.FullPairCountBudget(core, S, pair_polynomials, caps,
                                               basis_intervals, gauge, cuts)
    budget.establish_nonnegativity(result['pairs'])
    budget = refined_budget_helper.CancellationPreservingPairBudget(
        budget, verified_pair_majorants, result['pairs'])
    result['complete_pair_bound_method'] = budget.contract
    result['complete_pair_majorant_piece_count'] = len(verified_pair_majorants)
    result['full_pair_gauge_nonnegativity'] = 'ESTABLISHED_BY_THIS_RUNS_FRESH_EXACT_COMPLETE_BASE_PAIR_PROOFS'
    result['full_pair_gauge_actual_pair_multiplier'] = '0'
    result['full_pair_gauge_origin_constant'] = '0'
    result['full_pair_gauge_pair_identity'] = 'w=-q/(1-psi); q includes every original pair contribution; q+w*(1-psi)=0.'
    checkpoint()

    if stage in ['all', 'cyclic']:
        # No numerical status or pre-existing report enables any worker.
        # All these exact parent gates were executed in this invocation.
        assert result['parent_scalar_and_dependencies'] == 'FRESH_EXACT_REPLAY_PASSED'
        assert result['origin_identity'] == 'EXACT_MU_EQUALS_V_AND_VU_EQUALS_KAPPA'
        assert result['pointed_origin_cancellation'] == 'EXACT_CAP_TENT_CONDITIONAL_AND_TWO_CAPACITY_ANTIPODAL_LINEAR_ORIGIN_IDENTITIES'
        assert budget.nonnegative_proven
        assert_sources_unchanged()
        result['parallel']['parent_gates_before_fork'] = 'FRESH_DEPENDENCIES_SCALAR_ORIGIN_OBJECTIVE_AND_ALL_BASE_PAIRS_PASSED'
        result['identities'] = dict(status='INCOMPLETE', polynomial_identity_checks=0,
                                    angular_degree=angular_degree, branch_count=0)
        checkpoint()
        _WORKER_STATE = dict(
            parent_ready=True, folder=folder, core=core, budget=budget,
            identities=identities, angular_degree=angular_degree,
            expected_jobs=frozenset(jobs), cone=cone, coefficients=coefficients,
            flags=flags, joint_polynomials=joint_polynomials, conditional=conditional,
            tents=tents, basis_intervals=basis_intervals, lift_probes=lift_probes,
            caps=caps, lift_helper=lift_helper, allocations=allocations,
            spherical_rays=spherical_rays, graph_weights=graph_weights)
        expected, seen = set(jobs), set()
        executor = None
        futures = {}
        try:
            # Explicit fork is required: loaded modules, exact matrices and
            # the freshly proved budget are inherited, never serialized.
            context = multiprocessing.get_context('fork')
            executor = concurrent.futures.ProcessPoolExecutor(
                max_workers=workers, mp_context=context)
            futures = {executor.submit(_branch_worker, job): job for job in jobs}
            for future in concurrent.futures.as_completed(futures):
                job = futures[future]
                packet = future.result()
                result['parallel']['branch_packets'].append(packet)
                result['cyclic'].append(packet['cyclic'])
                # Save a returned VIOLATION/UNKNOWN before the acceptance
                # gate rejects it. Such a packet is never counted as proved.
                checkpoint()
                validate_branch_packet(job, packet, expected, seen, angular_degree)
                result['parallel']['accepted_jobs'] = len(seen)
                result['identities']['polynomial_identity_checks'] = 6*len(seen)
                result['identities']['branch_count'] = len(seen)
                print('FAST_PARALLEL_BRANCH_PASSED', job,
                      'completed', len(seen), 'of', len(jobs), flush=True)
                checkpoint()
            assert seen == expected
            identity_count = assert_complete_branch_coverage(
                jobs, result['parallel']['branch_packets'], angular_degree)
            order = {job: i for i, job in enumerate(jobs)}
            result['parallel']['branch_packets'].sort(
                key=lambda row: order[(row['kind'], tuple(row['regions']))])
            result['cyclic'].sort(
                key=lambda row: order[(row['kind'], tuple(row['regions']))])
            result['identities'].update(status='VERIFIED',
                                        polynomial_identity_checks=identity_count,
                                        branch_count=len(jobs))
            result['parallel']['coverage_status'] = 'ALL_EXPECTED_UNIQUE_BRANCHES_AND_IDENTITIES_VERIFIED'
        except BaseException as error:
            result['status'] = 'FAILED_OR_INCOMPLETE_EXACT_WIDEST_COORDINATE_BINARY_COMBINED_SIGNED_SPHERICAL_COUNT_GAUGE_FAST_PARALLEL_VERIFICATION'
            result['parallel']['failure'] = dict(type=type(error).__name__, message=str(error))
            checkpoint()
            # Stop pending exact jobs on a failed branch or interrupted run.
            # The processes are our own fork workers, never external jobs.
            for future in futures:
                future.cancel()
            if executor is not None:
                for process in list((getattr(executor, '_processes', None) or {}).values()):
                    if process.is_alive():
                        process.terminate()
            raise
        finally:
            try:
                if executor is not None:
                    executor.shutdown(wait=True, cancel_futures=True)
            finally:
                _WORKER_STATE = None
        checkpoint()

    # Repeat the arithmetic and complete source contract after every worker
    # has exited; no progress file or saved input status is an acceptance gate.
    assert objective == F(data['objective']) < F(data['threshold'])
    assert_sources_unchanged()
    assert split_backend.source_contract() == split_contract
    assert graph_objective == objective+graph_origin_constant/8 < F(graph_data['threshold'])
    result['final_input_sha256'] = {name: hashlib.sha256((folder/name).read_bytes()).hexdigest()
                                   for name in input_names+code_names}
    result['final_contract_sha256'] = {path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
                                      for path in contract_hashes}
    if stage != 'all':
        result['status'] = 'EXACT_WIDEST_COORDINATE_BINARY_COMBINED_SIGNED_SPHERICAL_COUNT_GAUGE_FAST_PARALLEL_'+stage.upper()+'_PROBE_PASSED_NOT_A_GLOBAL_CERTIFICATE'
        checkpoint()
        return result
    assert result['parallel']['coverage_status'] == 'ALL_EXPECTED_UNIQUE_BRANCHES_AND_IDENTITIES_VERIFIED'
    assert_complete_branch_coverage(jobs, result['parallel']['branch_packets'], angular_degree)
    assert result['identities']['status'] == 'VERIFIED'
    assert result['identities']['polynomial_identity_checks'] == 6*len(jobs)
    assert len(result['cyclic']) == len(jobs)
    result['status'] = 'EXACT_GLOBAL_BINARY_COMBINED_SIGNED_SPHERICAL_GRAPH_COUNT_GAUGE_CERTIFICATE_VERIFIED'
    result['input_sha256'] = {name: input_hashes[name] for name in input_names}
    result['verifier_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['seconds'] = time.time()-started
    _atomic_json(folder/'VERIFIED.json', result)
    checkpoint()
    print('ALL EXACT CONDITIONS PASSED: Delta_4 <', float(F(graph_data['threshold'])), flush=True)
    return result


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError('Run without -O.')
    parser = argparse.ArgumentParser()
    parser.add_argument('folder', type=Path)
    parser.add_argument('--stage', choices=['all', 'pair', 'cyclic'], default='all')
    parser.add_argument('--workers', type=int, choices=range(1, 5), default=1)
    args = parser.parse_args()
    verify(args.folder, args.stage, args.workers)
