"""Saved-certificate acceptance; no discovery, aliases, or historical branch trust.

Run from the project root. The unchanged global signs are a dependency on the
locally completed full V7 replay. All new local hypotheses are checked afresh.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse
import hashlib
import importlib.util
import json
import sys

PROJECT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
ROOT = HERE / 'source'
SOURCE = ROOT / 'current/candidates/round9_six_graphs_degree18'
OC = ROOT / 'target625/origin_completion'
GEO = ROOT / 'target625/geometry'
OLD = PROJECT / 'work/d4/2026-10-19-d4-v7-evidence-intake'
sys.dont_write_bytecode = True
sys.set_int_max_str_digits(0)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def read(path):
    def unique(pairs):
        out = {}
        for key, value in pairs:
            require(key not in out, 'duplicate JSON key')
            out[key] = value
        return out
    return json.loads(Path(path).read_text(), object_pairs_hook=unique)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def save(path, data):
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(data, indent=2) + '\n')
    temp.replace(path)


def stable(data):
    if isinstance(data, dict):
        return {k: stable(v) for k, v in data.items() if k != 'seconds'}
    if isinstance(data, list):
        return [stable(v) for v in data]
    return data


def check(output):
    require(__debug__ and not sys.flags.optimize, 'Python -O is forbidden')
    output.mkdir(parents=True, exist_ok=False)
    bindings = {}
    manifest = read(ROOT / 'MANIFEST.json')

    def bind(path, expected=None):
        path = Path(path).resolve()
        require(path.is_file() and path.is_relative_to(PROJECT), 'missing/escaping input')
        digest = sha(path)
        if path.is_relative_to(ROOT):
            row = manifest[path.relative_to(ROOT).as_posix()]
            require(path.stat().st_size == row['bytes'] and digest == row['sha256'], 'manifest mismatch')
        require(expected is None or digest == expected, 'binding mismatch: ' + str(path))
        bindings[str(path.relative_to(PROJECT))] = digest
        return path

    # Source identity is checked before importing any new checker.
    for path in ROOT.rglob('*.py'):
        bind(path)
    bind(Path(__file__))
    gate = load('local_v7_gate', bind(OLD / 'branches/accept_full_replay.py'))
    inherited = gate.check_run(OLD / 'rebuild/full', OLD / 'rebuild/full_receipt.json',
        source_map=read(OLD / 'branches/source_relocation.example.json'))
    require(inherited['cyclic_branches'] == 112, 'incomplete inherited signs')
    save(output / 'INHERITED_FULL_V7_ACCEPTANCE.json', inherited)
    accepted = read(bind(ROOT / 'audit/round9_six_graphs_monotone_origin_credit_global_verified.json'))
    # This old record supplies identity/accounting only, never branch acceptance.
    for name, digest in accepted['proof_file_sha256'].items():
        bind(ROOT / name, digest)
    verified_root = OLD / 'rebuild/full/record'
    for name in ('kernel.json', 'origin.json', 'scalar.json', 'graph_energies.json'):
        prior = bind(verified_root / 'current/candidates/round9_six_graphs_degree18' / name)
        bind(SOURCE / name, sha(prior))
    require(Q(accepted['objective']) == Q(inherited['objective']), 'old objective mismatch')
    origin = read(SOURCE / 'origin.json')
    kernel = read(SOURCE / 'kernel.json')
    graph_cost = Q(4769, 200000000)
    raw = Q(origin['objective']) + graph_cost
    old_credit = 8 * (raw - Q(inherited['objective']))
    require(raw == Q(accepted['exact_source_objective']) and old_credit > 0, 'raw accounting')

    projection_path = bind(OC / 'ACCEPTED_JOINT_COUNT_PROJECTION.json')
    projection = read(projection_path)
    projection_checker = load('local_projection', OC / 'verify_joint_count_origin_credit.py')
    fresh_projection = projection_checker.verify(projection, SOURCE)
    save(output / 'PROJECTION.json', fresh_projection)
    counts = load('local_counts', OC / 'fresh_joint_count_dependencies.py')
    # Only path resolution changes. This regular-file tree is the exact bound
    # historical subtree; the checker still verifies its original digests.
    counts.BASE = ROOT / 'research/sphere_packing_4d_research_record/historical_snapshot/new_work'
    fresh_counts = counts.verify(SOURCE,
        ROOT / 'audit/round9_six_graphs_monotone_origin_credit_global_verified.json', output / 'counts')
    require(fresh_counts['proved_count_bounds'] == [24, 24, 24, 25, 26], 'count caps')
    require(fresh_counts['proved_implications'] ==
        ['(23-N0)*(N4-N0)>=0', '(23-N1)*(N3-N1)>=0'], 'count implications')
    save(output / 'COUNT_GEOMETRY.json', fresh_counts)
    print('PASS projection and geometric count domain', flush=True)

    composer = load('reviewed_registry_helpers', OC / 'compose_radial_stratified_origin_credit.py')
    for path, digest in composer.PINS.items():
        bind(path, digest)
    registry = read(bind(OC / 'FINAL_ORIGIN_REGISTRY.json'))
    require(registry['schema'] == 'ACCEPTED_ORIGIN_RADIAL_AND_COMMON_G_ENVELOPE_V1', 'registry schema')
    for path, digest in registry['source_sha256'].items():
        bind(ROOT / path, digest)
    radial = load('local_radial', OC / 'verify_ray_center_radial_floor.py')
    radial_reports = []
    for entry in registry['radial_bundles']:
        bundle = read(bind(ROOT / entry['file'], entry['sha256']))
        require(len(bundle['exact_inputs']) == len(bundle['records']), 'radial record count')
        for exact, prior in zip(bundle['exact_inputs'], bundle['records']):
            fresh = radial.verify(exact, SOURCE, node_limit=100000)
            require(fresh == prior, 'radial replay differs')
            radial_reports.append(fresh)
    require(len(radial_reports) == 12, 'expected twelve radial planes')
    save(output / 'RADIAL_PLANES.json', radial_reports)
    rays = origin['appended_spherical_local_energy_rays']
    require(len(rays) == 2, 'spherical energy count')
    scale = sum(Q(ray['weight']) for ray in rays)
    B = sum(Q(ray['weight']) * Q(ray['B']) for ray in rays)
    g = [sum(Q(ray['weight']) * Q(ray['scalar_coefficients'][j]) for ray in rays) for j in range(13)]
    radii = list(map(Q, kernel['breaks'][1:6]))
    for ray in rays:
        require(Q(ray['weight']) > 0 and Q(ray['R']) == radii[-1] and
            Q(ray['cosine_cutoff']) == Q(21, 40) and Q(ray['global_origin']) == 0 and
            Q(ray['global_origin_coupling_coefficient']) == Q(ray['B'])/2 and
            list(map(Q, ray['global_kernel_angular_coefficients'])) == [-Q(v) for v in ray['scalar_coefficients']],
            'spherical energy physical identity')
    profile = [dict(N=N, spherical_surplus='0') for N in range(27)]
    stratified = []
    for index, entry in enumerate(registry['common_g_certificates']):
        path = bind(ROOT / entry['certificate'], entry['certificate_sha256'])
        prior = read(bind(ROOT / entry['prior_verification'], entry['prior_verification_sha256']))
        verifier_path, status = composer.KINDS[entry['kind']]
        module = load('new_spherical_' + str(index), verifier_path)
        print('CHECK angular', index, entry['kind'], flush=True)
        fresh = module.verify(path, node_limit=1200000)
        require(fresh['status'] == status and fresh['pair']['status'] == fresh['triangle']['status'] == 'CERTIFIED', 'angular failure')
        require(stable(fresh) == stable(prior), 'angular replay differs')
        combination = fresh['accepted_combination']
        require(Q(combination['positive_normalization']) == scale and
            Q(combination['weighted_accepted_B']) == B and
            list(map(Q, combination['weighted_accepted_g'])) == g and
            list(map(Q, fresh['common_scalar_coefficients'])) == [v/scale for v in g] and
            combination['accepted_origin_sha256'] == sha(SOURCE/'origin.json'), 'same-polynomial identity')
        if entry['kind'] == 'bounded_N26':
            require(fresh['maximum_code_count'] == 26 and Q(fresh['cosine_cutoff']) == Q(21,40), 'bounded scope')
            per = B - scale * Q(fresh['normalized_B'])
            require(per == Q(fresh['weighted_surplus_per_neighbor']) > 0, 'bounded surplus')
            for row in profile:
                row['spherical_surplus'] = str(max(Q(row['spherical_surplus']), row['N'] * per))
        else:
            require(entry['kind'] == 'radius_fixed_N', 'unexpected registry kind')
            N, j = fresh['fixed_count'], fresh['condition_radius_index']
            R, s = Q(fresh['condition_radial_radius']), Q(fresh['cosine_cutoff'])
            require(type(N) is int and N in (23,24,25) and type(j) is int and 0 <= j < 5, 'radius/count type')
            require(R == radii[j] and R*R < 2 and 1-Q(1)/(2*R*R) <= s <= Q(21,40), 'radius geometry')
            surplus = N * (B - scale * Q(fresh['B']))
            require(surplus == Q(fresh['weighted_surplus_at_N']) > 0, 'stratified surplus')
            stratified.append(dict(fixed_count=N, radius_index=j, surplus=str(surplus)))
        save(output / ('ANGULAR_' + str(index) + '.json'), fresh)
    # The latest bounded proof dominates all three older unstratified proofs.
    # We do not trust or rerun them: omission is safe. Equality of the final
    # credit below checks that the smaller, fully fresh proof suffices.
    enumerator = load('local_enumerator', OC / 'enumerate_radial_origin_envelope.py')
    local = enumerator.enumerate_profile(projection, radial_reports, profile, stratified)
    require(local['count_vectors_checked'] == 165503, 'incomplete enumeration')
    credit = Q(local['exact_origin_credit'])
    final = read(bind(ROOT / 'consolidation/FINAL_RESULT.json'))
    require(credit == Q(final['new_raw_origin_credit']) > old_credit, 'credit mismatch')
    objective = raw - credit/8
    threshold = Q(125907805731, 200000000000)
    require(objective == Q(final['exact_objective']) < threshold < Q(inherited['strict_upper_bound']), 'final objective')
    require(Q((objective*10**12).numerator//(objective*10**12).denominator+1, 10**12) == threshold, 'outward rounding')
    save(output / 'LOCAL_ENUMERATION.json', local)
    for name, digest in bindings.items():
        require(sha(PROJECT/name) == digest, 'input changed during replay')
    report = dict(status='LOCAL_FULL_V7_AND_FRESH_STRATIFIED_ORIGIN_ACCEPTED', dimension=4,
        packing_family='all sphere packings, including arbitrary upper density',
        global_density_bound_proved=True, exact_source_objective=str(raw),
        replaced_origin_credit=str(old_credit), new_origin_credit=str(credit),
        objective=str(objective), strict_upper_bound=str(threshold),
        prior_strict_upper_bound=inherited['strict_upper_bound'],
        threshold_improvement=str(Q(inherited['strict_upper_bound'])-threshold),
        count_vectors=165503, radial_planes=12, angular_proofs_fresh=3,
        complete_inherited_cyclic_branches=112, new_global_branch_run=False,
        source_bindings=bindings,
        trust_scope='Written local-energy proof; exact Python/NumPy integer arithmetic and outward Bernstein subdivision; existing full V7 replay dependency. No Lean or independent peer review claimed.')
    save(output / 'ACCEPTANCE.json', report)
    print('ACCEPTED strict all-packings D4 ceiling', threshold, flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    check(args.output.resolve())
