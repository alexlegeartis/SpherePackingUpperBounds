#!/usr/bin/env python3
"""Global origin completion from a registry of freshly verified local floors.

The inherited global candidate is fixed.  Registry data may supply exact
radial supporting planes and additional proofs for its SAME weighted
spherical energy.  All proofs are checked by pinned verifier families;
their cardinality and radius hypotheses are enforced during enumeration.
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
GEO = ROOT/'target625/geometry'
SOURCE = ROOT/'current/candidates/round9_six_graphs_degree18'
INHERITED = HERE/'compose_uniform_joint_count_global_credit.py'
RADIAL = HERE/'verify_ray_center_radial_floor.py'
ENUMERATOR = HERE/'enumerate_radial_origin_envelope.py'
PROJECTION = HERE/'ACCEPTED_JOINT_COUNT_PROJECTION.json'
FIXED = GEO/'certify_accepted_combined_fixedN_credit.py'
BOUNDED = GEO/'certify_accepted_combined_boundedN_credit.py'
STRATIFIED = GEO/'certify_accepted_stratified_fixedN_credit.py'
CORE = GEO/'certify_fixedN_common_spherical.py'
FIXED_ALGEBRA = GEO/'fixed_cardinality_spherical_algebra.py'
FULL_ALGEBRA = ROOT/'optimization/continuation_20261018/base/new_work/spherical_sdp/certify.py'
PINS = {
    INHERITED: 'dc7fe375db41a81bde8535739862d997c8b9487bb676634594fe1826200eb9a3',
    RADIAL: '0528ced8fbcdd52966c77de700d0b2de4200745f7959a7b40d825da628eafd5c',
    ENUMERATOR: '52c10655c95f37bbb540cd095184b037e2b4a5ca34798a8681432f0efec52906',
    PROJECTION: '219da1d0c222b941552b42899528bbdc4a86074d2fb49ab2ea1aeb7d13cdaa7c',
    FIXED: '664e13edbef85e0664b1d6732f32d8be800561026be5d1e2cd76e3df82f02a6c',
    BOUNDED: '9296ccbfd4468de4f2738b78b78d6625eced56e5ef5a2d4f46498c9e81e75cfb',
    STRATIFIED: '338732eb88bd0af68dc3b3d5a9bc2e3b200da896dfe6786a3cde55d9162249fc',
    CORE: '03e5584dfb1d19740035342fa7afb0a497b8283e6c8ae4ddf27303a5e2c21536',
    FIXED_ALGEBRA: '66174183e7ef69ebf496035e05b9b16e1875faefdfe45cec25566efa3997d5ed',
    FULL_ALGEBRA: 'aeb873edff74112cd0dd0293d1c77be91aec7ec405aa1013b9b92088db46fca5',
}
KINDS = {
    'fixed_N': (FIXED, 'EXACT_FIXED_CARDINALITY_SURPLUS_OF_ACCEPTED_SPHERICAL_SUM_VERIFIED'),
    'bounded_N26': (BOUNDED, 'EXACT_BOUNDED_CARDINALITY_SURPLUS_OF_ACCEPTED_SPHERICAL_SUM_VERIFIED'),
    'radius_fixed_N': (STRATIFIED, 'EXACT_RADIUS_STRATIFIED_FIXED_CARDINALITY_SURPLUS_OF_ACCEPTED_SPHERICAL_SUM_VERIFIED'),
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
        return {k: stable(v) for k, v in value.items() if k != 'seconds'}
    if isinstance(value, list):
        return [stable(v) for v in value]
    return value


def relative_file(name):
    path = Path(name)
    assert not path.is_absolute() and '..' not in path.parts
    result = (ROOT/path).resolve()
    assert result.is_relative_to(ROOT.resolve()) and result.is_file()
    return result


def compose(registry_path, output_dir):
    assert __debug__, 'Run without -O.'
    started = time.time()
    registry_path = Path(registry_path).resolve()
    assert registry_path.is_relative_to(ROOT.resolve())
    output_dir = Path(output_dir).resolve()
    assert output_dir.is_relative_to(ROOT.resolve())
    output_dir.mkdir(parents=True, exist_ok=False)
    binding = {str(path.relative_to(ROOT)): digest for path, digest in PINS.items()}

    def bind(path, digest=None):
        path = Path(path)
        actual = sha(path)
        if digest is not None:
            assert actual == digest, path
        key = str(path.relative_to(ROOT))
        assert key not in binding or binding[key] == actual, path
        binding[key] = actual
        return path

    for path, digest in PINS.items():
        assert sha(path) == digest, path
    own = Path(__file__)
    bind(own)
    bind(registry_path)
    registry = json.loads(registry_path.read_text())
    assert registry['schema'] == 'ACCEPTED_ORIGIN_RADIAL_AND_COMMON_G_ENVELOPE_V1'
    expected_sources = {str((SOURCE/name).relative_to(ROOT)): sha(SOURCE/name)
                        for name in ('kernel.json', 'origin.json')}
    assert registry['source_sha256'] == expected_sources
    for path, digest in expected_sources.items():
        bind(ROOT/path, digest)
    radial_entries = registry['radial_bundles']
    spherical_entries = registry['common_g_certificates']
    assert isinstance(radial_entries, list) and isinstance(spherical_entries, list)
    assert len({row['file'] for row in radial_entries}) == len(radial_entries)
    assert len({row['certificate'] for row in spherical_entries}) == len(spherical_entries)
    for row in radial_entries:
        bind(relative_file(row['file']), row['sha256'])
    for row in spherical_entries:
        assert row['kind'] in KINDS
        bind(relative_file(row['certificate']), row['certificate_sha256'])
        bind(relative_file(row['prior_verification']), row['prior_verification_sha256'])
    inherited_module = load('radial_registry_inherited_composer', INHERITED)
    inherited, inherited_path = inherited_module.compose(output_dir/'inherited_uniform_replay')
    assert inherited['status'] == 'EXACT_GLOBAL_FIVE_COUNT_AND_MAX_COMMON_G_ORIGIN_CREDIT_VERIFIED'
    assert inherited['global_density_bound_proved'] is True
    assert inherited['original142_payloads_unchanged'] is True
    assert inherited['complete_domain_coverage']['total_branches'] == 112
    assert inherited['count_geometry_freshly_replayed'] is True
    assert inherited['boundedN26_proof_freshly_replayed'] is True
    assert inherited['fixedN_proofs_freshly_replayed'] == 2
    assert inherited['all_pair_and_cyclic_expressions_unchanged'] is True
    for path, digest in inherited['proof_file_sha256'].items():
        bind(ROOT/path, digest)
    bind(inherited_path)
    radial_module = load('radial_registry_exact_plane_checker', RADIAL)
    radial_reports = []
    for index, row in enumerate(radial_entries):
        bundle = json.loads(relative_file(row['file']).read_text())
        assert bundle['status'] == 'EXACT_RAY_CENTER_RADIAL_ORIGIN_FLOOR_BUNDLE_VERIFIED'
        assert bundle['global_density_bound_proved'] is False
        assert bundle['verifier_sha256'] == PINS[RADIAL]
        assert len(bundle['exact_inputs']) == len(bundle['records']) > 0
        reports = []
        for exact, prior in zip(bundle['exact_inputs'], bundle['records']):
            fresh = radial_module.verify(exact, SOURCE, node_limit=100000)
            assert fresh == prior
            assert fresh['all_seven_closed_radial_pieces_verified'] is True
            reports.append(fresh)
            for path, digest in fresh['source_sha256'].items():
                bind(ROOT/path, digest)
        radial_reports += reports
        fresh_path = output_dir/f'FRESH_RADIAL_BUNDLE_{index}.json'
        write(fresh_path, dict(status='FRESH_EXACT_RADIAL_ORIGIN_FLOORS_VERIFIED', records=reports))
        bind(fresh_path)
    origin = json.loads((SOURCE/'origin.json').read_text())
    kernel = json.loads((SOURCE/'kernel.json').read_text())
    radii = list(map(Q, kernel['breaks'][1:6]))
    rays = origin['appended_spherical_local_energy_rays']
    assert len(rays) == 2
    scale = sum(Q(ray['weight']) for ray in rays)
    B = sum(Q(ray['weight'])*Q(ray['B']) for ray in rays)
    g = [sum(Q(ray['weight'])*Q(ray['scalar_coefficients'][j]) for ray in rays) for j in range(13)]
    assert scale > 0
    for ray in rays:
        assert Q(ray['weight']) > 0 and Q(ray['R']) == radii[-1]
        assert Q(ray['cosine_cutoff']) == Q(21, 40)
        assert Q(ray['global_origin']) == 0
        assert Q(ray['global_origin_coupling_coefficient']) == Q(ray['B'])/2
        assert list(map(Q, ray['global_kernel_angular_coefficients'])) == [-Q(x) for x in ray['scalar_coefficients']]
    spherical_profile = [dict(N=row['N'], spherical_surplus=row['spherical_surplus'])
                         for row in inherited['conditional_count_profile']]
    assert [row['N'] for row in spherical_profile] == list(range(27))
    stratified = []
    spherical_reports = []
    for index, row in enumerate(spherical_entries):
        path = relative_file(row['certificate'])
        verifier, expected_status = KINDS[row['kind']]
        module = load(f'radial_registry_common_g_{index}', verifier)
        if row['kind'] == 'bounded_N26':
            assert Path(module.algebra.__file__).resolve() == FULL_ALGEBRA.resolve()
            assert sha(Path(module.algebra.__file__)) == PINS[FULL_ALGEBRA]
        else:
            assert sha(Path(module.core.__file__)) == PINS[CORE]
            assert sha(Path(module.core.algebra.__file__)) == PINS[FIXED_ALGEBRA]
        if hasattr(module, 'common_record'):
            assert sha(Path(module.common_record.__code__.co_filename)) == PINS[FIXED]
        fresh = module.verify(path, node_limit=1200000)
        assert fresh['status'] == expected_status
        assert fresh['global_bound_proved'] is False
        assert fresh['pair']['status'] == fresh['triangle']['status'] == 'CERTIFIED'
        assert fresh['certificate_sha256'] == row['certificate_sha256']
        assert stable(fresh) == stable(json.loads(relative_file(row['prior_verification']).read_text()))
        combination = fresh['accepted_combination']
        assert combination['accepted_origin_sha256'] == sha(SOURCE/'origin.json')
        assert Q(combination['positive_normalization']) == scale
        assert Q(combination['weighted_accepted_B']) == B
        assert list(map(Q, combination['weighted_accepted_g'])) == g
        assert list(map(Q, fresh['common_scalar_coefficients'])) == [value/scale for value in g]
        if row['kind'] == 'bounded_N26':
            assert fresh['maximum_code_count'] == 26
            assert fresh['common_g_unchanged'] is True and Q(fresh['cosine_cutoff']) == Q(21, 40)
            surplus = B-scale*Q(fresh['normalized_B'])
            assert surplus == Q(fresh['weighted_surplus_per_neighbor']) > 0
            for item in spherical_profile:
                item['spherical_surplus'] = str(max(Q(item['spherical_surplus']), item['N']*surplus))
        else:
            N = fresh['fixed_count']
            assert type(N) is int and N in (23, 24, 25)
            surplus = N*(B-scale*Q(fresh['B']))
            assert surplus == Q(fresh['weighted_surplus_at_N']) > 0
            if row['kind'] == 'fixed_N':
                assert Q(fresh['cosine_cutoff']) == Q(21, 40)
                spherical_profile[N]['spherical_surplus'] = str(max(Q(spherical_profile[N]['spherical_surplus']), surplus))
            else:
                j = fresh['condition_radius_index']
                R = Q(fresh['condition_radial_radius'])
                s = Q(fresh['proof_cosine_cutoff'])
                assert type(j) is int and 0 <= j < 5 and radii[j] == R
                # The derivative numerator r^2-s^2+1 is >=2-R^2>0;
                # hence both radii maximize their direction cosine at R.
                assert 1 <= R <= radii[-1] and 2-R*R > 0
                rho = 1-Q(1)/(2*R*R)
                assert Q(fresh['geometric_cosine_upper']) == rho <= s <= Q(21, 40)
                assert Q(fresh['cosine_cutoff']) == s
                stratified.append(dict(fixed_count=N, radius_index=j, surplus=str(surplus),
                                       condition='N_j=N equals the total fifth near count, so every supported point has radius<=R_j.'))
        spherical_reports.append(fresh)
        fresh_path = output_dir/f'FRESH_COMMON_G_SURPLUS_{index}.json'
        write(fresh_path, fresh)
        bind(fresh_path)
    projection = json.loads(PROJECTION.read_text())
    # The inherited pinned composer has freshly checked this exact data.
    assert projection['count_vectors_checked'] == 165503
    enumerator = load('radial_registry_count_enumerator', ENUMERATOR)
    local = enumerator.enumerate_profile(projection, radial_reports, spherical_profile, stratified)
    assert local['status'] == 'EXACT_RADIAL_PROJECTION_MAXIMUM_COUNT_ENUMERATION_COMPLETE'
    assert local['count_vectors_checked'] == 165503
    assert [row['N'] for row in local['conditional_count_profile']] == list(range(27))
    credit = Q(local['exact_origin_credit'])
    assert credit == min(Q(row['exact_floor']) for row in local['conditional_count_profile'])
    previous = Q(inherited['new_origin_credit'])
    original_credit = Q(inherited['old_origin_credit'])
    assert credit > previous > original_credit > 0
    raw = Q(inherited['exact_source_objective'])
    graph = Q(inherited['graph_origin_constant_retained'])
    assert raw == Q(origin['objective'])+graph/8 == Q(inherited['objective'])+previous/8
    objective = raw-credit/8
    scaled = objective*10**12
    threshold = Q(scaled.numerator//scaled.denominator+1, 10**12)
    assert objective < threshold < Q(inherited['strict_upper_bound'])
    local.update(all_external_dependencies_freshly_verified=True,
                 exact_radial_records_freshly_verified=len(radial_reports),
                 additional_common_g_certificates_freshly_verified=len(spherical_reports),
                 stratified_spherical_surpluses=stratified,
                 unstratified_spherical_surplus_profile=spherical_profile,
                 aggregate_local_energy='Original E0 + the three original count energies + the two original spherical energies - this new total credit >=0.',
                 origin_credit_accounting='Use max of E0 bounds and max of spherical bounds, add the count energies once, and replace every older origin credit once.')
    local_path = output_dir/'LOCAL_RADIAL_STRATIFIED_ORIGIN_CREDIT_VERIFIED.json'
    write(local_path, local)
    bind(local_path)
    for path, digest in binding.items():
        assert sha(ROOT/path) == digest, path
    record = dict(status='EXACT_GLOBAL_RADIAL_FEASIBILITY_AND_STRATIFIED_COMMON_G_ORIGIN_CREDIT_VERIFIED',
                  global_density_bound_proved=True, source=str(SOURCE.relative_to(ROOT)),
                  registry_file=str(registry_path.relative_to(ROOT)), registry_sha256=sha(registry_path),
                  inherited_fresh_global_proof=str(inherited_path.relative_to(ROOT)),
                  inherited_fresh_global_proof_sha256=sha(inherited_path),
                  original_accepted_global_proof_sha256=inherited['original_accepted_global_proof_sha256'],
                  original142_payloads_unchanged=True,
                  complete_domain_coverage=inherited['complete_domain_coverage'],
                  all_pair_and_cyclic_expressions_unchanged=True,
                  all_graphs_retained_unchanged=True, nonorigin_kernel_unchanged=True,
                  total_physical_origin_coupling_unchanged=True,
                  graph_origin_constant_retained=str(graph),
                  graph_density_cost_retained=inherited['graph_density_cost_retained'],
                  old_origin_credit=str(original_credit), replaced_previous_credit=str(previous),
                  new_origin_credit=str(credit), exact_density_credit=str(credit/8),
                  additional_density_improvement=str((credit-original_credit)/8),
                  improvement_over_previous_density=str((credit-previous)/8),
                  exact_source_objective=str(raw), objective=str(objective), objective_decimal=float(objective),
                  strict_upper_bound=str(threshold), threshold_decimal=float(threshold),
                  count_vectors_checked=165503, conditional_count_profile=local['conditional_count_profile'],
                  radial_origin_floors_freshly_verified=len(radial_reports),
                  additional_common_g_certificates_freshly_verified=len(spherical_reports),
                  local_certificate=str(local_path.relative_to(ROOT)), local_certificate_sha256=sha(local_path),
                  origin_credit_accounting=local['origin_credit_accounting'],
                  theorem='Exact ray-square completions and radial polynomial bounds strengthen the E0 lower bound. Source-matched common-g theorems strengthen the same spherical energy only under their verified count/radius hypotheses. Complete integer enumeration gives a uniform local credit; the unchanged complete global sign certificate and Fourier normalization then imply the stated strict density bound for arbitrary four-dimensional sphere packings.',
                  new_global_interval_run_required=False, new_global_interval_run_claimed=False,
                  proof_file_sha256=binding, composer_source=str(own.relative_to(ROOT)), composer_sha256=sha(own),
                  elapsed_seconds=time.time()-started)
    output = output_dir/'GLOBAL_RADIAL_STRATIFIED_ORIGIN_CREDIT_VERIFIED.json'
    write(output, record)
    return record, output


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--registry', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    record, output = compose(args.registry, args.output_dir)
    print(json.dumps(dict(status=record['status'], strict_upper_bound=record['strict_upper_bound'],
                         objective_decimal=record['objective_decimal'], output=str(output),
                         output_sha256=sha(output)), indent=2))
