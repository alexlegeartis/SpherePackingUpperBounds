#!/usr/bin/env python3
"""Independent kernel check: compare kernel.py with stage 1 exactly.

    python3 run.py --out DIR [--workers N]

DIR must not exist. The run
  1. checks this folder against MANIFEST.json and the stage 1 record and
     data against stage 1's MANIFEST.json;
  2. copies verify/stage1_v7/record to DIR/work/record and places the four
     data/v7 files where stage 1's run.py places them;
  3. runs extract_stage1.py there: the unchanged stage 1 coordinator runs its
     parent gates and fresh pair proofs, and every one of the 112 cyclic
     cases is compared exactly with kernel.py (see README.md);
  4. checks A, U_raw, the 7 pair pieces and the 112 cases, and runs the
     internal consistency checks of kernel.py;
  5. writes DIR/RESULT.json. Exit code 0 if and only if status is PASS.

Negative-control options (they must give FAIL):
  --perturb implementation   changes one spherical weight inside kernel.py
  --kernel-data DIR          kernel.py reads this data folder instead of
                             data/v7 (stage 1 still reads data/v7)
Testing option: --only 'full:0,0,0;outside:5,5' compares only these cases
(the result is then FAIL for incompleteness).
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import itertools
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
STAGE1 = PKG/'verify'/'stage1_v7'
CENSUS = PKG/'data'/'census_floors.json'


class Fail(Exception):
    pass


def require(cond, message):
    if not cond:
        raise Fail(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_own_manifest():
    manifest = json.loads((HERE/'MANIFEST.json').read_text())
    for row in manifest['files']:
        path = HERE/row['path']
        require(path.is_file() and sha(path) == row['sha256'],
                'kernel_independent file differs from MANIFEST.json: '+row['path'])


def check_stage1_manifest():
    manifest = json.loads((STAGE1/'MANIFEST.json').read_text())
    for row in manifest['files']:
        if not row['path'].startswith('record/'):
            continue
        require(sha(STAGE1/row['path']) == row['sha256'], 'stage 1 file differs: '+row['path'])
    for row in manifest['data_files']:
        require(sha((STAGE1/row['path']).resolve()) == row['sha256'],
                'data file differs from stage 1 MANIFEST.json: '+row['path'])
    return manifest


def self_consistency(K):
    """kernel.py's polynomial form equals its own pointwise definition.

    For every case, the homogenized polynomial is evaluated at two interior
    rational points and compared with factor * (sum of the polynomial
    components of K at the three anchors), evaluated pointwise from the
    Bernstein formula, Chebyshev recurrence and indicator definitions.
    """
    import flint
    checks = 0
    for kind, size in [('full', 3), ('outside', 2)]:
        for regions in itertools.combinations_with_replacement(range(K.np), size):
            poly = K.gamma_poly(kind, regions)
            for frac in [(Q(1, 3), Q(2, 5), Q(5, 7)), (Q(6, 7), Q(1, 9), Q(1, 2))]:
                point = [K.breaks[p]+(K.breaks[p+1]-K.breaks[p])*f for p, f in zip(regions, frac)]
                if kind == 'outside':
                    point.append(frac[2])
                a, b, c = point
                fv = poly(*(flint.fmpq(x.numerator, x.denominator) for x in (a, b, c)))
                value = Q(int(fv.p), int(fv.q))
                if kind == 'full':
                    cos = lambda u, v, w: (u*u+v*v-w*w)/(2*u*v)
                    direct = (K.pair_kernel_point(a, b, cos(a, b, c))+K.pair_kernel_point(a, c, cos(a, c, b))
                              + K.pair_kernel_point(b, c, cos(b, c, a)))
                else:
                    direct = K.pair_kernel_point(a, b, K.outside_cosine(a, b, c))
                require(value == K.factor_point(kind, point)*direct,
                        'kernel.py self-consistency failed at %s %s' % (kind, regions))
                checks += 1
    # pair: D(r) + truncation + cap self term equals 2h + k(x,x) pointwise
    for p in range(K.np):
        D = K.pair_poly(p)
        lo, hi = K.breaks[p:p+2]
        for f in [Q(1, 3), Q(4, 5)]:
            r = lo+(hi-lo)*f
            value = sum((c*r**a for a, c in enumerate(D)), Q(0))
            direct = 2*K.h_point(r)+K.pair_kernel_point(r, r, 1)+K.cap_self_diag_point(r)
            # D uses untruncated cap linear terms; J_equivalent removes them exactly.
            require(value-K.J_equivalent_point(r) == direct, 'pair self-consistency failed on piece %d' % p)
            checks += 1
    return checks


def origin_identity(K):
    """Gram origin block: M_0 u = v and u.v = kappa >= 0 (49 equations)."""
    n = K.n
    for i in range(n):
        row = sum((Q(int(K.M[0][i, j]), K.D)*K.u[j] for j in range(n)), Q(0))
        require(row == K.v[i], 'M u != v at %d' % i)
    kappa = sum((a*b for a, b in zip(K.u, K.v)), Q(0))
    require(kappa == K.kappa and kappa >= 0, 'u.v != kappa')
    return n


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--workers', type=int, choices=[1, 2], default=1)
    parser.add_argument('--only', default=None)
    parser.add_argument('--perturb', choices=['implementation'], default=None)
    parser.add_argument('--kernel-data', type=Path, default=None)
    parser.add_argument('--cxx', default=os.environ.get('CXX', 'clang++'))
    args = parser.parse_args()
    if not __debug__:
        parser.error('run without -O')
    out = args.out.expanduser().resolve()
    if out.exists():
        parser.error('--out must be a new directory')
    out.mkdir(parents=True)
    started = time.time()
    result = dict(stage='kernel_independent', status='FAIL')
    kernel_data = (args.kernel_data or DATA).resolve()
    result['kernel_data'] = str(kernel_data)
    result['negative_control'] = dict(perturb=args.perturb,
                                      kernel_data_is_package_data=kernel_data == DATA.resolve())
    try:
        check_own_manifest()
        manifest = check_stage1_manifest()
        sys.path.insert(0, str(HERE))
        import kernel as kmod
        import flint
        import numpy
        result['runtime'] = dict(python=platform.python_version(), python_flint=flint.__version__,
                                 numpy=numpy.__version__, platform=platform.platform())
        K = kmod.Kernel(kernel_data)
        if args.perturb == 'implementation':
            K.sph[0]['weight'] *= 1+Q(1, 10**9)
        work = out/'work'
        rec = work/'record'
        shutil.copytree(STAGE1/'record', rec, ignore=shutil.ignore_patterns('__pycache__'))
        for row in manifest['data_files']:
            target = rec/row['placed_at']
            require(not target.exists(), 'data file already present: '+row['placed_at'])
            shutil.copyfile((STAGE1/row['path']).resolve(), target)
            require(sha(target) == row['sha256'], 'copied data file differs')
        cxx = shutil.which(args.cxx)
        require(cxx is not None, 'C++ compiler not found: '+args.cxx)
        toolchain = out/'toolchain'
        toolchain.mkdir()
        (toolchain/'g++').write_text('#!/bin/sh\nexec "'+cxx+'" "$@"\n')
        (toolchain/'g++').chmod(0o755)
        env = {k: v for k, v in os.environ.items() if k != 'PYTHONOPTIMIZE'}
        env.update(PYTHONDONTWRITEBYTECODE='1', OPENBLAS_NUM_THREADS='1',
                   PATH=str(toolchain)+os.pathsep+os.environ.get('PATH', ''))
        cases = out/'cases'
        argv = [sys.executable, str(HERE/'extract_stage1.py'), '--record', str(rec),
                '--data', str(kernel_data), '--kernel-dir', str(HERE), '--out', str(cases),
                '--workers', str(args.workers)]
        if args.only:
            argv += ['--only', args.only]
        if args.perturb:
            argv += ['--perturb', args.perturb]
        t0 = time.time()
        with (out/'extract.log').open('w') as log:
            code = subprocess.run(argv, cwd=rec, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
        result['extraction'] = dict(exit=code, seconds=round(time.time()-t0, 1))
        require(code == 0, 'extract_stage1.py exited %d; see extract.log' % code)

        # A and U_raw.
        U = K.U_raw()
        A = K.A()
        coordinator = json.loads((cases/'coordinator.json').read_text())
        result['coordinator'] = coordinator
        U1 = Q(coordinator['stage1_objective'])
        A1 = 8*U1-K.F[0]
        census = Q(json.loads(CENSUS.read_text())['exact_raw_objective'])
        result.update(A=str(A), U_raw=str(U), U_raw_decimal='%.15f' % float(U),
                      A_parts={k: str(v) for k, v in K.origin_parts().items()},
                      A_equal_stage1=(A == A1), U_raw_equal_stage1=(U == U1),
                      U_raw_equal_census=(U == census))
        failures = []
        contracts_path = cases/'contracts.json'
        result['component_contracts'] = (json.loads(contracts_path.read_text())
                                        if contracts_path.is_file() else None)
        if not result['component_contracts'] or result['component_contracts'].get('status') != 'PASS':
            failures.append('component parameter and boundary contracts missing or failed')
        if not (A == A1 and U == U1 and U == census):
            failures.append('A or U_raw differs from stage 1 or census_floors.json')
        if coordinator.get('stage1_pairs_certified') != 7:
            failures.append('stage 1 did not freshly certify its 7 pair pieces in this run')
        if coordinator.get('stage1_origin_identity') != 'EXACT_MU_EQUALS_V_AND_VU_EQUALS_KAPPA':
            failures.append('stage 1 origin identity gate missing')

        # Pair pieces.
        pair = json.loads((cases/'pair.json').read_text()) if (cases/'pair.json').is_file() else None
        result['pair'] = pair
        if pair is None or pair['status'] != 'PASS' or len(pair['pieces']) != 7:
            failures.append('pair comparison missing or failed')

        # Cyclic cases.
        expected = [(kind, regions) for kind, size in [('full', 3), ('outside', 2)]
                    for regions in itertools.combinations_with_replacement(range(K.np), size)]
        rows, totals = [], dict(full=0, outside=0, terms=0, sample_points=0)
        outcomes = {}
        for kind, regions in expected:
            path = cases/('%s_%s.json' % (kind, '_'.join(map(str, regions))))
            if not path.is_file():
                rows.append(dict(kind=kind, regions=list(regions), status='MISSING'))
                continue
            row = json.loads(path.read_text())
            rows.append(dict(kind=kind, regions=list(regions), status=row['status'],
                             terms=row.get('polynomial', {}).get('stage1_terms'),
                             stage1_sha256=row.get('polynomial', {}).get('stage1_sha256'),
                             eta=row.get('eta'), seconds=row.get('seconds'),
                             failure=(row.get('error') or row.get('polynomial', {}).get('mismatches')
                                      or row.get('samples', {}).get('failures') or None)))
            if row['status'] == 'PASS':
                totals[kind] += 1
                totals['terms'] += row['polynomial']['stage1_terms']
                totals['sample_points'] += row['samples']['points']
                for comp, how in row['samples']['outcomes'].items():
                    for k, v in how.items():
                        outcomes.setdefault(comp, {}).setdefault(k, 0)
                        outcomes[comp][k] += v
        result['cases'] = rows
        result['counts'] = dict(cases_expected=len(expected), full_equal=totals['full'],
                                outside_equal=totals['outside'],
                                polynomial_terms_compared=totals['terms'],
                                pair_pieces_equal=sum(r['D_equal'] and r['F_equal'] and not r['J_mismatches']
                                                      for r in (pair or {}).get('pieces', [])),
                                correction_sample_points=totals['sample_points'],
                                correction_sample_outcomes=outcomes)
        if totals['full'] != 84 or totals['outside'] != 28 or len(expected) != 112:
            bad = [r for r in rows if r['status'] != 'PASS']
            failures.append('%d of 112 cases not equal; first: %s' % (len(bad), json.dumps(bad[:2])[:2000]))

        # kernel.py internal checks.
        result['self_consistency_checks'] = self_consistency(K)
        result['origin_identity_equations'] = origin_identity(K)
        result['identity_scope'] = {
            'polynomial_parts': 'Exact coefficientwise identity on all 112 closed branch polynomials',
            'pair_F_D': 'Exact univariate coefficients on all 7 pieces',
            'origin_objective': 'Exact rational identity',
            'component_parameters': 'Exact equality of validated inputs and 49 radial basis polynomials',
            'nonpolynomial_formulas': 'Written proof in COMPONENT_IDENTITIES.md; samples are diagnostics only',
            'pair_J_samples': 'Diagnostic samples on rational subintervals; not a standalone continuum proof',
        }
        result['failures'] = failures
        if not failures:
            result['status'] = 'PASS'
    except Exception as error:
        result['status'] = 'FAIL'
        result['failure'] = '%s: %s' % (type(error).__name__, error)
    result['seconds'] = round(time.time()-started, 1)
    (out/'RESULT.json').write_text(json.dumps(result, indent=1)+'\n')
    print(json.dumps({k: result.get(k) for k in ('stage', 'status', 'counts', 'U_raw_decimal',
                                                 'failure', 'failures', 'seconds')}, indent=1)[:4000])
    return 0 if result['status'] == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
