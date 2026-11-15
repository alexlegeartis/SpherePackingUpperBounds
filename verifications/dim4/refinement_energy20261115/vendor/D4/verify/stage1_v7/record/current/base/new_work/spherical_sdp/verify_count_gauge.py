"""Exact global verification with counts, pointed balls, tents and exclusions.

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

No numerical solver output, failed branch or incomplete stage is acceptance.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import importlib.util
import itertools
import json
import math
import sys
import time
import numpy as np
from count_quadratic import NestedCountCone

HERE = Path(__file__).resolve().parent


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


def verify(folder, stage='all'):
    if not __debug__:
        raise RuntimeError('Run without -O.')
    assert stage in ['all', 'pair', 'cyclic']
    folder = Path(folder).resolve()
    started = time.time()
    input_names = ['origin.json', 'scalar.json', 'kernel.json']
    code_names = ['certificate_core.py', 'verify_scalar.py', 'audit_identities.py']
    input_hashes = {name: hashlib.sha256((folder/name).read_bytes()).hexdigest()
                    for name in input_names+code_names}
    contract_paths = [Path(__file__), HERE/'count_quadratic.py', HERE/'certify.py',
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
    sys.path.insert(0, str(folder))
    core = load_module('certificate_core', folder/'certificate_core.py')
    scalar = load_module('verify_scalar', folder/'verify_scalar.py')
    assert core.ROOT.resolve() == folder
    data = json.loads((folder/'origin.json').read_text())
    S = json.loads((folder/'scalar.json').read_text())
    # The scalar proof must cover the entire complement of kernel support.
    # A larger scalar-tail cutoff would leave an unverified radial gap.
    assert 1 <= F(S['cutoff']) <= core.BREAKS[-1]
    kernel_data = json.loads((folder/'kernel.json').read_text())
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
    assert 'local_energy' not in data and 'spherical_local_energy' not in data
    sys.path.insert(0, str(HERE.parent/'pointed'))
    pointed_path = HERE.parent/'pointed'/'interval_certificate_full_closed.py'
    pointed = load_module('pointed_general_full_exact', pointed_path)
    pair_pointed = load_module('pointed_general_pair_exact', HERE.parent/'pointed'/'pair_interval_certificate_closed.py')
    tent_helper = load_module('geometric_tent_exact', HERE.parent/'pointed'/'interval_certificate_tent_closed.py')
    budget_helper = load_module('full_pair_count_gauge_exact', HERE.parent/'pointed'/'full_pair_count_gauge.py')
    lift_helper = load_module('count_gauge_cyclic_exact', HERE.parent/'pointed'/'interval_certificate_count_gauge.py')
    affine = load_module('geometric_affine_exact', HERE.parent/'conditional'/'affine_latitude_counts.py')
    joint_helper = load_module('geometric_joint_count_exact', HERE.parent/'conditional'/'joint_count_polynomial.py')
    assert data['pointed_caps'], 'This contract requires explicit pointed corrections.'
    basis_intervals = [tuple(map(F, row))
                       for row in kernel_data['radial_basis_intervals']]
    caps = pointed.validate_caps(data['pointed_caps'], core, basis_intervals)
    tents = tent_helper.validate_tents(data.get('pointed_tents', []), core, basis_intervals)
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
    result = {'status': 'INCOMPLETE_EXACT_FULL_PAIR_COUNT_GAUGE_VERIFICATION',
              'count_polynomial': count_polynomial, 'pairs': [], 'cyclic': [],
              'pointed_cap_count': len(caps),
              'pointed_tent_count': len(tents),
              'full_pair_count_gauge': {key: [str(v) for v in value] if isinstance(value, list) else str(value)
                                        for key, value in gauge.items()},
              'tent_capacity_lemmas': [tent['capacity_lemma'] for tent in tents],
              'conditional_rays': conditional_reports,
              'conditional_count_polynomials': joint_reports,
              'truncated_pointed_cap_count': sum('pivot_lower_cutoff' in cap for cap in caps),
              'pointed_contract': 'GENERAL_SIGMA_AND_ACTIVE_NONORIGIN_BALL_COUNTS',
              'pointed_helper_sha256': hashlib.sha256(pointed_path.read_bytes()).hexdigest(),
              'initial_input_sha256': input_hashes,
              'contract_sha256': contract_hashes}
    progress_path = folder/('COUNT_GAUGE_PROGRESS.json' if stage == 'all' else
                            'COUNT_GAUGE_'+stage.upper()+'_PROBE.json')

    def checkpoint():
        result['seconds'] = time.time()-started
        progress_path.write_text(json.dumps(result, indent=2)+'\n')

    checkpoint()
    if stage == 'all':
        spherical = load_module('spherical_nested_exact', HERE/'certify.py')
        checked, proof_results = {}, []
        # Each copy is byte-identical; replay outputs cannot overwrite older
        # accepted reports.  The originals are hashed again at completion.
        replay_folder = folder/'spherical_replay'
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
        scalar.verify()
        if 'radial_basis_intervals' in kernel_data:
            identities = load_module('audit_identities', folder/'audit_identities.py')
            result['identities'] = identities.verify(folder)
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
        region = i//core.NC
        origin_shift[i] += sum((F(23, 2)*ray['weight']*conditional_flags(core, region, ray)[1]
                               for ray in conditional), F(0))
    assert all(v[i]+origin_shift[i] == pair_v[i] for i in range(n)), \
        'The auxiliary pair row must equal the PSD row plus every exact geometric linear coefficient.'
    result['pointed_origin_cancellation'] = 'EXACT_CAP_TENT_AND_CONDITIONAL_LINEAR_ORIGIN_IDENTITIES'
    pair_polynomials = {}
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
                record = pair_pointed.check(core, pc, pd, (lo, hi), region, caps,
                                             basis_intervals, alpha)
                record.update(region=region, left=str(lo), right=str(hi),
                              nested_pair_correction=str(correction),
                              conditional_polynomial_pair_correction=str(joint_pair),
                              count_flags=flags(region), active_linear_weight=str(active),
                              exponential_order=order, positive_part_bound=str(positive_bound))
                result['pairs'].append(record)
                print('GEOMETRIC_POINTED_PAIR', record, flush=True)
                checkpoint()
                assert record['status'] == 'CERTIFIED', record

    # Even a cyclic-only probe must freshly prove every complete unlifted
    # pair branch. A status claimed in an input file cannot enable the lift.
    budget = budget_helper.FullPairCountBudget(core, S, pair_polynomials, caps,
                                               basis_intervals, gauge, cuts)
    budget.establish_nonnegativity(result['pairs'])
    result['full_pair_gauge_nonnegativity'] = 'ESTABLISHED_BY_THIS_RUNS_FRESH_EXACT_COMPLETE_BASE_PAIR_PROOFS'
    result['full_pair_gauge_actual_pair_multiplier'] = '0'
    result['full_pair_gauge_origin_constant'] = '0'
    result['full_pair_gauge_pair_identity'] = 'w=-q/(1-psi); q includes every original pair contribution; q+w*(1-psi)=0.'
    checkpoint()

    if stage in ['all', 'cyclic']:
        for kind, number in [('full', 3), ('outside', 2)]:
            for regions in itertools.combinations_with_replacement(range(core.NS), number):
                active_flags = [flags(region) for region in regions]
                if kind == 'outside':
                    active_flags.append([0]*cone.k)
                correction = dot(cone.triangle_row_from_flags(active_flags), coefficients)
                joint_cyclic = F(0)
                for polynomial, weight in joint_polynomials:
                    joint_flags = [polynomial.flags(core.BREAKS, region) for region in regions]
                    if kind == 'outside':
                        joint_flags.append([0]*polynomial.k)
                    joint_cyclic += weight*dot(polynomial.triangle_row(joint_flags), polynomial.coefficients)
                correction += joint_cyclic
                mixed = conditional_correction(core, kind, regions, conditional)
                eta = -correction-mixed
                record = (tent_helper.check(core, kind, regions, eta, tents,
                                            basis_intervals, caps=caps) if budget is None else
                          lift_helper.check(core, kind, regions, eta, tents,
                                            basis_intervals, budget, lift_probes, caps=caps))
                record.update(kind=kind, eta=str(eta),
                              pointed_correction='EXACT_CAP_TENT_AND_CONDITIONAL_BOUNDS',
                              conditional_polynomial_cyclic_correction=str(joint_cyclic),
                              conditional_cyclic_correction=str(mixed))
                record.update(nested_cyclic_correction=str(correction),
                              count_flags=active_flags)
                result['cyclic'].append(record)
                print('GEOMETRIC_POINTED_CYCLIC', record, flush=True)
                checkpoint()
                assert record['status'] == 'CERTIFIED', record

    objective = (fs[0]+kappa+coefficients[0]+
                 sum((cut['M']*F(cut['weight']) for cut in cuts), F(0))+
                 sum((weight*polynomial.coefficients[0]
                      for polynomial, weight in joint_polynomials), F(0)))/8
    assert alpha == 2 and F(S['G'][0]) == 1
    assert objective == F(data['objective']) < F(data['threshold'])
    assert all(hashlib.sha256((folder/name).read_bytes()).hexdigest() == value
               for name, value in input_hashes.items()), \
        'Certificate inputs changed during verification.'
    assert all(hashlib.sha256(Path(path).read_bytes()).hexdigest() == value
               for path, value in contract_hashes.items()), \
        'Verifier sources changed during verification.'
    assert all(hashlib.sha256((HERE/proof['certificate']).read_bytes()).hexdigest()
               == proof['certificate_sha256'] for proof in all_dependencies), \
        'A spherical count certificate changed during verification.'
    result.update(objective=str(objective), objective_decimal=float(objective),
                  strict_upper_bound=data['threshold'])
    if stage != 'all':
        result['status'] = 'EXACT_COUNT_GAUGE_'+stage.upper()+'_PROBE_PASSED_NOT_A_GLOBAL_CERTIFICATE'
        checkpoint()
        return result
    result['status'] = 'EXACT_GLOBAL_FULL_PAIR_COUNT_GAUGE_CERTIFICATE_VERIFIED'
    result['input_sha256'] = {name: input_hashes[name] for name in input_names}
    result['verifier_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['seconds'] = time.time()-started
    (folder/'VERIFIED.json').write_text(json.dumps(result, indent=2)+'\n')
    checkpoint()
    print('ALL EXACT CONDITIONS PASSED: Delta_4 <', float(F(data['threshold'])), flush=True)
    return result


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError('Run without -O.')
    parser = argparse.ArgumentParser()
    parser.add_argument('folder', type=Path)
    parser.add_argument('--stage', choices=['all', 'pair', 'cyclic'], default='all')
    args = parser.parse_args()
    verify(args.folder, args.stage)
