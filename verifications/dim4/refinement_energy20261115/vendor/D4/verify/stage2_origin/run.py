"""Stage 2: radial and angular origin credit and the 165503-vector count census.

Public driver adapted from the private acceptance driver verify.py. Only path
plumbing and the stage 1 hand-off differ; see CHANGES.md. Every local proof is
recomputed from exact inputs. Stored files are used only as inputs or as
expected values that the fresh results must equal.

Usage:
    python3 run.py --out DIR [--stage1 STAGE1_DIR]

DIR must not exist. Exit 0 if and only if DIR/RESULT.json has status PASS.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse
import hashlib
import importlib.util
import json
import shutil
import sys
import time
import traceback

sys.dont_write_bytecode = True
sys.set_int_max_str_digits(0)

STAGE = Path(__file__).resolve().parent            # verify/stage2_origin
PKG = STAGE.parents[1]                             # the D4 folder
V7 = PKG / 'data/v7'                               # kernel, origin, scalar, graph_energies
EXPORT = STAGE / 'source'                          # exported subset of the private source tree
REFERENCE = STAGE / 'reference'                    # accepted output of the private run (2026-10-20)
CANDIDATE = 'current/candidates/round9_six_graphs_degree18'
V7_NAMES = ('kernel.json', 'origin.json', 'scalar.json', 'graph_energies.json')
# Entries of the old accepted record's proof_file_sha256 that this stage reads.
# The other entries are inputs of the stage 1 V7 sign proof (see CHANGES.md).
HIST = 'research/sphere_packing_4d_research_record/historical_snapshot/new_work'
STAGE2_PROOF_FILES = [CANDIDATE + '/' + name for name in V7_NAMES] + [
    HIST + '/spherical_sdp/certify.py',
    HIST + '/spherical_sdp/count_quadratic.py',
    HIST + '/spherical_sdp/s05075_taper12_exact.json',
    HIST + '/spherical_sdp/s0515_count25_exact.json',
    HIST + '/spherical_sdp/s0525_count26_exact.json',
    HIST + '/conditional/affine_latitude_counts.py',
    HIST + '/conditional/joint_count_polynomial.py',
]

# Set by build_tree(): the run tree is a fresh copy, so that every original
# module resolves its ROOT = Path(__file__).parents[2] inside it.
PROJECT = ROOT = SOURCE = OC = None


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


def build_tree(out):
    """Copy the exported source and the four V7 data files into OUT/tree."""
    global PROJECT, ROOT, SOURCE, OC
    PROJECT = out / 'tree'
    ROOT = PROJECT / 'source'
    SOURCE = ROOT / CANDIDATE
    OC = ROOT / 'target625/origin_completion'
    shutil.copytree(EXPORT, ROOT, ignore=shutil.ignore_patterns('__pycache__', '.DS_Store'))
    SOURCE.mkdir(parents=True)
    for name in V7_NAMES:
        shutil.copyfile(V7 / name, SOURCE / name)


def stage1_inherited(stage1):
    """Replacement for the private gate.check_run on the local full V7 replay.

    With --stage1, the stage 1 output must report PASS, its exact raw objective
    must equal the one recomputed here, its record must hold the same four V7
    files, and its old origin composition supplies the inherited objective.
    Without --stage1, the committed stage 1 reference values are used.
    """
    reference = read(REFERENCE / 'INHERITED_FULL_V7_ACCEPTANCE.json')
    if stage1 is None:
        return dict(objective=reference['objective'], strict_upper_bound=reference['strict_upper_bound'],
                    cyclic_branches=reference['cyclic_branches'], stage1_checked=False,
                    source='reference/INHERITED_FULL_V7_ACCEPTANCE.json (committed stage 1 values)')
    result = read(stage1 / 'RESULT.json')
    require(result.get('stage') == 'stage1_v7' and result.get('status') == 'PASS',
            'stage 1 RESULT.json status is not PASS')
    require(result.get('mode') == 'full' and result.get('branches_rerun') == 112,
            'stage 1 did not rerun all 112 cyclic branches (quick mode is not enough)')
    record = stage1 / 'record'
    for name in V7_NAMES:
        require(sha(record / CANDIDATE / name) == sha(SOURCE / name), 'stage 1 record differs from data/v7: ' + name)
    composition = read(record / 'replay/v7_full_origin_composition.json')
    require(composition['status'] == 'EXACT_GLOBAL_ORIGIN_CREDIT_COMPOSITION_VERIFIED'
            and composition['global_density_bound_proved'] is True
            and composition['complete_global_route'] ==
            'NORMAL_V7_ALL_112_BRANCHES; NO_HISTORICAL_PARTIAL_PROOF_COMPOSITION', 'stage 1 composition incomplete')
    require(Q(composition['objective']) == Q(result['credited_objective'])
            and Q(composition['strict_upper_bound']) == Q(result['strict_upper_bound']),
            'stage 1 RESULT.json and composition differ')
    return dict(objective=composition['objective'], strict_upper_bound=composition['strict_upper_bound'],
                cyclic_branches=result['branches_rerun'], stage1_checked=True,
                stage1_raw_objective=result['raw_objective'], source=str(stage1))


def check(output, stage1):
    require(__debug__ and not sys.flags.optimize, 'Python -O is forbidden')
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
    inherited = stage1_inherited(stage1)
    require(inherited['cyclic_branches'] == 112, 'incomplete inherited signs')
    save(output / 'STAGE1_INHERITED.json', inherited)
    accepted = read(bind(ROOT / 'audit/round9_six_graphs_monotone_origin_credit_global_verified.json'))
    # This old record supplies identity/accounting only, never branch acceptance.
    for name in STAGE2_PROOF_FILES:
        bind(ROOT / name, accepted['proof_file_sha256'][name])
    require(Q(accepted['objective']) == Q(inherited['objective']), 'old objective mismatch')
    origin = read(SOURCE / 'origin.json')
    kernel = read(SOURCE / 'kernel.json')
    graph_cost = Q(4769, 200000000)
    raw = Q(origin['objective']) + graph_cost
    old_credit = 8 * (raw - Q(inherited['objective']))
    require(raw == Q(accepted['exact_source_objective']) and old_credit > 0, 'raw accounting')
    if inherited['stage1_checked']:
        require(Q(inherited['stage1_raw_objective']) == raw, 'stage 1 raw objective differs')

    projection_path = bind(OC / 'ACCEPTED_JOINT_COUNT_PROJECTION.json')
    projection = read(projection_path)
    projection_checker = load('local_projection', OC / 'verify_joint_count_origin_credit.py')
    fresh_projection = projection_checker.verify(projection, SOURCE)
    save(output / 'PROJECTION.json', fresh_projection)
    counts = load('local_counts', OC / 'fresh_joint_count_dependencies.py')
    # Only path resolution changes. This regular-file tree is the exact bound
    # historical subtree; the checker still verifies its original digests.
    counts.BASE = ROOT / HIST
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
    return report, inherited


# Fresh outputs must equal the committed reference except for timings
# ('seconds'). In ACCEPTANCE.json the source_bindings keys differ by prefix
# only: 'source/X' here, 'work/d4/radial-stratified-origin-bound/source/X' in
# the private run. Every fresh binding must name the same bytes as the
# reference binding of X.
PRIVATE_PREFIX = 'work/d4/radial-stratified-origin-bound/'


def comparable(data):
    if isinstance(data, dict):
        return {k: comparable(v) for k, v in data.items() if k not in ('seconds', 'source_bindings')}
    if isinstance(data, list):
        return [comparable(v) for v in data]
    return data


def compare_with_reference(out):
    names = ['ANGULAR_0.json', 'ANGULAR_1.json', 'ANGULAR_2.json', 'COUNT_GEOMETRY.json',
             'PROJECTION.json', 'RADIAL_PLANES.json', 'LOCAL_ENUMERATION.json', 'ACCEPTANCE.json']
    names += ['counts/' + p.name for p in sorted((REFERENCE / 'counts').iterdir()) if p.suffix == '.json']
    differences = [name for name in names if comparable(read(out / name)) != comparable(read(REFERENCE / name))]
    fresh = read(out / 'ACCEPTANCE.json')['source_bindings']
    prior = read(REFERENCE / 'ACCEPTANCE.json')['source_bindings']
    require(fresh and all(prior.get(PRIVATE_PREFIX + key) == digest for key, digest in fresh.items()),
            'a bound input differs from the input of the accepted private run')
    return names, differences, len(fresh)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--out', type=Path, required=True, help='new output directory')
    parser.add_argument('--stage1', type=Path, help='output directory of verify/stage1_v7/run.py')
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    started = time.time()
    result = dict(stage='stage2_origin', status='FAIL')
    try:
        stage1 = args.stage1.resolve() if args.stage1 else None
        build_tree(out)
        report, inherited = check(out, stage1)
        names, differences, bound = compare_with_reference(out)
        require(not differences, 'fresh output differs from reference: ' + ', '.join(differences))
        enumeration = read(out / 'LOCAL_ENUMERATION.json')
        groups = enumeration['first_containing_radius_profile']
        require(len(groups) == 124, 'expected 124 (N, first radius) groups')
        result.update(
            status='PASS',
            stage1_checked=inherited['stage1_checked'], stage1_source=inherited['source'],
            exact_source_objective=report['exact_source_objective'],
            inherited_objective=inherited['objective'],
            replaced_origin_credit=report['replaced_origin_credit'],
            new_origin_credit=report['new_origin_credit'],
            objective=report['objective'], strict_upper_bound=report['strict_upper_bound'],
            count_vectors_checked=enumeration['count_vectors_checked'], group_minima=len(groups),
            radial_planes=report['radial_planes'], angular_proofs_fresh=report['angular_proofs_fresh'],
            spherical_count_bounds=[24, 24, 24, 25, 26],
            bound_inputs_identical_to_accepted_run=bound,
            matches_reference=names)
    except Exception as error:  # fail closed on any error
        traceback.print_exc()
        result['error'] = '%s: %s' % (type(error).__name__, error)
    result['seconds'] = round(time.time() - started, 1)
    save(out / 'RESULT.json', result)
    print(json.dumps({k: v for k, v in result.items() if k != 'matches_reference'}, indent=2))
    return 0 if result['status'] == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
