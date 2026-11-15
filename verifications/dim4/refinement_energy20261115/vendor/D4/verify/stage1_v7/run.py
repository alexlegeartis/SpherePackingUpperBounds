#!/usr/bin/env python3
"""Stage 1 of the D4 bound: rerun the V7 base certificate from exported sources.

    python3 run.py --out DIR [--workers N] [--quick]

Full mode (default) recomputes every exact check, including all 112 cyclic
triangle branches, then runs the full origin composer. Quick mode runs the
same parent checks and pair pieces, the spherical count and scalar checks,
and the composer in --prepare mode, but it does NOT rerun the 112 cyclic
branches. RESULT.json says which mode ran.

DIR must not exist. Layout of DIR, as `replay.py --mode full --output DIR`:
  DIR/record/...   working copy of the exported record, with the V7 data
                   files taken from ../../data/v7/
  DIR/replay.log   output of every step
  DIR/checks/      outputs of the independent checkers
  DIR/RESULT.json  final status
Exit code 0 if and only if RESULT.json has status PASS.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
PKG = HERE.parents[1]
DATA = PKG/'data'/'v7'
CENSUS = PKG/'data'/'census_floors.json'
CANDIDATE = 'current/candidates/round9_six_graphs_degree18'
COORDINATOR = 'current/kernel_integration_route/verify_combined_binary_six_exterior_graphs_small_cap_bounds.py'
COMPOSER = 'audit/compose_six_graph_full_origin_credit.py'
CREDIT = 'audit/round9_origin_credit/ORIGIN_CREDIT_VERIFIED.json'
GLOBAL_STATUS = 'EXACT_GLOBAL_BINARY_COMBINED_SIGNED_SPHERICAL_GRAPH_COUNT_GAUGE_CERTIFICATE_VERIFIED'
PAIR_STATUS = 'EXACT_WIDEST_COORDINATE_BINARY_COMBINED_SIGNED_SPHERICAL_COUNT_GAUGE_FAST_PARALLEL_PAIR_PROBE_PASSED_NOT_A_GLOBAL_CERTIFICATE'
PAIR_FILE = 'SIX_EXTERIOR_GRAPH_SMALL_CAP_BOUNDS_CORRELATED_FACTOR_WIDEST_COORDINATE_BINARY_COMBINED_SIGNED_SPHERICAL_COUNT_GAUGE_FAST_PARALLEL_PAIR_PROBE.json'
COMPOSED_STATUS = 'EXACT_GLOBAL_ORIGIN_CREDIT_COMPOSITION_VERIFIED'
PREPARED_STATUS = 'EXACT_ORIGIN_CREDIT_COMPOSITION_PREPARED_NOT_GLOBAL'
GRAPH_CHARGE = Q(4769, 200000000)
THRESHOLD = Q(78749932299, 125000000000)


class Fail(Exception):
    pass


def require(condition, message):
    if not condition:
        raise Fail(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def exact(value):
    require(type(value) in (str, int), 'exact value must be a string or integer')
    return Q(value)


def check_manifest():
    """Every exported file and every V7 data file must match the manifest bytes."""
    manifest = json.loads((HERE/'MANIFEST.json').read_text())
    for row in manifest['files']:
        path = HERE/row['path']
        require(path.is_file(), 'missing exported file: '+row['path'])
        data = path.read_bytes()
        require(len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256'],
                'exported file differs from MANIFEST.json: '+row['path'])
    listed = {row['path'] for row in manifest['files']}
    # Finder may add this metadata while a user browses the package. It is
    # never executable or a certificate input; all other unlisted files fail.
    present = {str(p.relative_to(HERE)) for p in (HERE/'record').rglob('*')
               if p.is_file() and p.name != '.DS_Store'}
    require(present <= listed, 'unlisted files in record/: '+', '.join(sorted(present-listed)[:5]))
    for row in manifest['data_files']:
        path = (HERE/row['path']).resolve()
        require(path.is_file(), 'missing data file: '+row['path'])
        require(sha(path) == row['sha256'], 'data file differs from MANIFEST.json: '+row['path'])
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--out', type=Path, required=True, help='new, nonexistent output directory')
    parser.add_argument('--workers', type=int, choices=range(1, 5), default=2,
                        help='parallel cyclic workers (the coordinator allows 1 to 4)')
    parser.add_argument('--quick', action='store_true', help='skip the 112 cyclic branches')
    parser.add_argument('--cxx', default=os.environ.get('CXX', 'clang++'),
                        help='C++20 compiler used for the native helpers (default clang++)')
    args = parser.parse_args()
    if not __debug__:
        parser.error('run without -O; the checkers rely on assert')
    out = args.out.expanduser().resolve()
    if out.exists():
        parser.error('--out must be a new directory')
    if out.is_relative_to(HERE):
        parser.error('--out must be outside the exported tree')
    mode = 'quick' if args.quick else 'full'
    started = time.time()
    out.mkdir(parents=True)
    log_path = out/'replay.log'
    steps = []
    result = {'stage': 'stage1_v7', 'status': 'FAIL', 'mode': mode,
              'branches_rerun': 0, 'workers': args.workers}

    def run(name, argv, cwd):
        print('RUN', name, flush=True)
        t0 = time.time()
        with log_path.open('a') as log:
            log.write('\n=== RUN '+name+': '+' '.join(map(str, argv))+'\n')
            log.flush()
            code = subprocess.run(argv, cwd=cwd, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
        steps.append({'name': name, 'argv': [str(a) for a in argv], 'exit': code,
                      'seconds': round(time.time()-t0, 1)})
        require(code == 0, f'step {name} exited {code}; see replay.log')

    try:
        import numpy
        import flint
        manifest = check_manifest()
        work = out/'record'
        shutil.copytree(HERE/'record', work, ignore=shutil.ignore_patterns('__pycache__'))
        for row in manifest['data_files']:
            target = work/row['placed_at']
            require(not target.exists(), 'data file already present in record: '+row['placed_at'])
            shutil.copyfile((HERE/row['path']).resolve(), target)
            require(sha(target) == row['sha256'], 'copied data file differs')
        (work/'replay').mkdir(exist_ok=True)
        checks = out/'checks'
        checks.mkdir()

        # The native helpers call 'g++'. Route that name to the chosen compiler.
        cxx = shutil.which(args.cxx)
        require(cxx is not None, 'C++ compiler not found: '+args.cxx)
        toolchain = out/'toolchain'
        toolchain.mkdir()
        shim = toolchain/'g++'
        shim.write_text('#!/bin/sh\nexec "'+cxx+'" "$@"\n')
        shim.chmod(0o755)
        compiler = subprocess.run([cxx, '--version'], capture_output=True, text=True).stdout.splitlines()
        env = {k: v for k, v in os.environ.items() if k != 'PYTHONOPTIMIZE'}
        env.update(PYTHONDONTWRITEBYTECODE='1', OPENBLAS_NUM_THREADS='1',
                   PATH=str(toolchain)+os.pathsep+os.environ.get('PATH', ''))
        result['runtime'] = {'python': platform.python_version(), 'numpy': numpy.__version__,
                             'python_flint': flint.__version__, 'flint_version': getattr(flint, '__FLINT_VERSION__', None),
                             'platform': platform.platform(), 'cxx': cxx, 'cxx_version': compiler[:1]}
        py = sys.executable

        # 1. The original V7 coordinator. Full: all stages, 112 branches. Quick: pair stage.
        stage = 'pair' if args.quick else 'all'
        run('coordinator --stage '+stage, [py, COORDINATOR, CANDIDATE, '--stage', stage,
                                           '--workers', str(args.workers)], work)
        folder = work/CANDIDATE
        origin = json.loads((folder/'origin.json').read_text())
        base = exact(origin['objective'])
        u_raw = base+GRAPH_CHARGE
        if args.quick:
            probe = json.loads((folder/PAIR_FILE).read_text())
            require(probe['status'] == PAIR_STATUS, 'pair stage status: '+str(probe['status']))
            report = probe
            run('quick parent checks', [py, HERE/'checks/quick_parent_checks.py', work,
                                        '--output', checks/'quick_parent_checks.json'], work)
            require(json.loads((checks/'quick_parent_checks.json').read_text())['status']
                    == 'QUICK_PARENT_CHECKS_PASSED', 'quick parent checks failed')
        else:
            verified = folder/'VERIFIED.json'
            require(verified.is_file(), 'coordinator wrote no VERIFIED.json')
            report = json.loads(verified.read_text())
            require(report['status'] == GLOBAL_STATUS, 'global status: '+str(report['status']))
            par = report['parallel']
            jobs = {(p['kind'], tuple(p['regions'])) for p in par['expected_jobs']}
            got = [(c['kind'], tuple(c['regions'])) for c in report['cyclic']]
            require(len(jobs) == 112 and len(got) == len(set(got)) == 112 and set(got) == jobs,
                    'branch coverage is not exactly the 112 expected branches')
            require(sum(k == 'full' for k, _ in jobs) == 84, 'expected 84 full branches')
            require(all(c['status'] == 'CERTIFIED' and not c['failures'] for c in report['cyclic']),
                    'a cyclic branch is not CERTIFIED')
            require(par['accepted_jobs'] == 112 and not par.get('unproved_jobs'), 'unaccepted branches')
            require(par['coverage_status'] == 'ALL_EXPECTED_UNIQUE_BRANCHES_AND_IDENTITIES_VERIFIED', 'coverage')
            require(report['identities']['polynomial_identity_checks'] == 672, 'identity count')
            require(report['parent_scalar_and_dependencies'] == 'FRESH_EXACT_REPLAY_PASSED', 'scalar gate')
            result['branches_rerun'] = 112
            result['cyclic'] = {'branches': 112, 'full': 84, 'exterior': 28,
                                'nodes': sum(int(c['nodes']) for c in report['cyclic']),
                                'leaves': sum(int(c['leaves']) for c in report['cyclic']),
                                'correlated_factor_leaves': par.get('correlated_factor_leaves'),
                                'polynomial_identity_checks': 672}
        # Gates shared by both modes: base pairs, origin identity, objective.
        pairs = report['pairs']
        require(len(pairs) == 7 and all(p['status'] == 'CERTIFIED' and not p['failures'] for p in pairs),
                'pair pieces not all CERTIFIED')
        require(report['origin_identity'] == 'EXACT_MU_EQUALS_V_AND_VU_EQUALS_KAPPA', 'origin identity')
        require(exact(report['base_objective']) == base, 'base objective differs from origin.json')
        require(exact(report['objective']) == u_raw, 'coordinator objective is not base + 4769/200000000')
        require(report['independent_graph_local_energies']['status'] == 'FRESH_EXACT_INDEPENDENT_GRAPH_ENERGIES_VERIFIED',
                'graph energies')
        result['pairs'] = {'pieces': len(pairs), 'nodes': sum(int(p['nodes']) for p in pairs),
                           'leaves': sum(int(p['leaves']) for p in pairs)}

        # 2. The original full origin composer (--prepare in quick mode).
        composition = 'replay/v7_full_origin_composition.json' if not args.quick else 'replay/v7_prepared_origin_composition.json'
        argv = [py, COMPOSER, CANDIDATE, CREDIT, composition]+(['--prepare'] if args.quick else [])
        run('composer'+(' --prepare' if args.quick else ''), argv, work)
        comp = json.loads((work/composition).read_text())
        if args.quick:
            require(comp['status'] == PREPARED_STATUS and comp['global_density_bound_proved'] is False,
                    'prepared composition status')
        else:
            require(comp['status'] == COMPOSED_STATUS and comp['global_density_bound_proved'] is True,
                    'composition status')
            require(comp['inherited_global_proof']['verified_sha256'] == sha(work/CANDIDATE/'VERIFIED.json'),
                    'composer did not bind this run\'s VERIFIED.json')
            require(exact(comp['strict_upper_bound']) == THRESHOLD, 'strict bound differs')
        require(exact(comp['source_objective']) == u_raw, 'composer source objective differs')
        credited = exact(comp['objective'])
        require(credited < THRESHOLD, 'credited objective misses threshold')

        # 3. Independent checkers written by reviewers (standard library only).
        run('check_v7_inputs', [py, HERE/'checks/check_v7_inputs.py', work,
                                '--output', checks/'v7_inputs.json'], work)
        require(json.loads((checks/'v7_inputs.json').read_text())['status']
                == 'V7_EXACT_DATA_GRAM_ORIGIN_AND_SCALAR_CHECKS_PASSED', 'v7 input check')
        run('audit_scalar_origin', [py, HERE/'checks/audit_scalar_origin.py', '--record', work,
                                    '--output', checks/'scalar_origin_supplement.json'], work)
        supplement = json.loads((checks/'scalar_origin_supplement.json').read_text())
        require(supplement['status'] == 'PASS_EXACT_PHYSICAL_ORIGIN_OBJECTIVE_STURM_AND_REJECTION_AUDIT'
                and supplement['negative_control_count'] == 11, 'scalar/origin supplement')
        rec = supplement['objective_reconstruction']
        require(exact(rec['base_objective'])+exact(rec['graph_charge']) == u_raw, 'supplement raw objective')
        require(exact(rec['final_objective']) == credited, 'supplement credited objective')
        run('verify_moderate_independent', [py, HERE/'checks/verify_moderate_independent.py',
            work/'current/six_exterior_graph_route/moderate_exterior_graph_certificate.json',
            '--out', checks/'moderate_independent.json'], work)
        moderate = json.loads((checks/'moderate_independent.json').read_text())
        require(moderate['status'] == 'INDEPENDENT_EXACT_MODERATE_GRAPH_CAPACITY_1880_PASS'
                and moderate['even_ordered_edge_capacity'] == 1880, 'moderate graph capacity')
        run('check_broad19_rays', [py, HERE/'checks/check_broad19_rays.py', work,
                                   '--output', checks/'broad19_rays.json'], work)
        require(json.loads((checks/'broad19_rays.json').read_text())['status']
                == 'BROAD19_RAYS_ARE_KERNEL_GRAM_RAYS', 'broad19 rays')
        run('replay_graphs', [py, HERE/'checks/replay_graphs.py', '--record', work,
                              '--output-dir', checks/'graphs'], work)
        graphs = json.loads((checks/'graphs'/'ACCEPTANCE.json').read_text())
        require(graphs['status'] == 'GRAPH_LOCAL_ACCEPTANCE_PASS'
                and graphs['malformed_input_rejections'] == 11
                and graphs['singleton_controls'] == 1406
                and graphs['rational_box_point_controls'] == 270
                and exact(graphs['combined_graph_density_charge']) == GRAPH_CHARGE, 'graph local acceptance')
        require([c.get('even_ordered_edge_bound') for c in graphs['commands'][:4]] == [1880, 120, 1358, 604],
                'graph capacities')

        # 4. Cross-check against the value used downstream by the census stage.
        if CENSUS.is_file():
            census = exact(json.loads(CENSUS.read_text())['exact_raw_objective'])
            require(census == u_raw, 'U_raw differs from data/census_floors.json exact_raw_objective')
            result['matches_census_exact_raw_objective'] = True

        result.update(
            status='PASS',
            raw_objective=str(u_raw), raw_objective_decimal=f'{float(u_raw):.15f}',
            raw_objective_definition='(F(0)+A)/8 = base objective + graph charge 4769/200000000, before the origin credit',
            base_objective=str(base), graph_charge=str(GRAPH_CHARGE),
            credited_objective=str(credited), credited_objective_decimal=f'{float(credited):.15f}',
            origin_credit_density=comp['density_credit'],
            strict_upper_bound=str(THRESHOLD) if not args.quick else None,
            global_sign_proof=('RECOMPUTED: all 112 cyclic branches and 7 pair pieces' if not args.quick else
                               'NOT RERUN: quick mode skips the 112 cyclic branches; pair pieces were recomputed'),
            checks={'v7_inputs': 'PASS', 'scalar_origin_supplement': 'PASS',
                    'moderate_independent': 'PASS', 'broad19_rays': 'PASS',
                    'graph_local_acceptance': 'PASS'} |({'quick_parent_checks': 'PASS'} if args.quick else {}))
    except Exception as error:
        result['status'] = 'FAIL'
        result['failure'] = f'{type(error).__name__}: {error}'
    result['steps'] = steps
    result['seconds'] = round(time.time()-started, 1)
    (out/'RESULT.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result.get(k) for k in ('stage', 'status', 'mode', 'branches_rerun',
                      'raw_objective_decimal', 'credited_objective_decimal', 'failure', 'seconds')
                      if k in result}, indent=2))
    return 0 if result['status'] == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
