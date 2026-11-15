"""Extract stage 1's exact kernel objects on the exterior branches over B_5.

Outer call (from run.py):  prepare(PKG, WORK) builds a working copy of the stage 1
record exactly as stage1_v7/run.py does (manifest-checked bytes, the four V7 data
files placed at their manifest paths, a g++ shim for the native helpers).

Inner call (subprocess, cwd = WORK/record):
  python3 -B extract_stage1.py --inner --out PICKLE

1. Runs the unmodified stage 1 coordinator verify(folder, 'pair', 1): fresh scalar,
   Gram, origin identity, graph capacities, count polynomials, the 7 base pair
   proofs and the full pair budget. A wrapper around the coordinator's load_module
   snapshots verify's local variables when it constructs the refined pair budget
   (its last object). Nothing else in verify is changed.
2. For every exterior branch ('outside', (i, j)), 0 <= i <= j <= 4, calls the
   coordinator's own _branch_worker with the captured objects. Only the final
   inequality check (lift_helper.check), the identity audit and the C++ factor
   backend are replaced by stubs that record their arguments. This yields the
   exact eta (count, conditional-polynomial and conditional-ray cyclic terms) and
   the exact correction objects (caps, tents, allocations, spherical rays, graph
   weights, budget) that the real check receives.
3. Calls the real build_signed(core, 'outside', regions, eta, spherical_rays), the
   first step of lift_helper.check: the integer numerator C with
       C(r, s, z) / dd = (2 r s)^18 * (Gram + signed spherical - eta),
   t = -1 + z((r + s)^2 - 49/16)/(2 r s).
4. Builds the subtracted corrections of stage 1 (gauge, caps, tents, allocations,
   graphs) on each ordered cell as exact symbolic values in (r, s, t), from stage 1's
   own formulas: the inner lambdas of radial_bounds, probability_bounds,
   outside_tent_bounds and squared_distance_bounds are called directly (taken from
   their code objects, closures rebuilt); the remaining formulas are transcribed
   with file and line given in the comments, and checked against stage 1's own
   bound routines at exact singleton points (singleton_checks).
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import importlib.util
import itertools
import json
import os
import pickle
import shutil
import subprocess
import sys
import time
import types

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
CANDIDATE = 'current/candidates/round9_six_graphs_degree18'
COORDINATOR = 'current/kernel_integration_route/verify_combined_binary_six_exterior_graphs_small_cap_bounds.py'
PAIR_STATUS = 'EXACT_WIDEST_COORDINATE_BINARY_COMBINED_SIGNED_SPHERICAL_COUNT_GAUGE_FAST_PARALLEL_PAIR_PROBE_PASSED_NOT_A_GLOBAL_CERTIFICATE'
NREG = 5  # B_5 uses radial regions 0..4


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ------------------------------------------------------------------ outer
def prepare(pkg, work, cxx='clang++'):
    stage1 = Path(pkg)/'verify'/'stage1_v7'
    manifest = json.loads((stage1/'MANIFEST.json').read_text())
    for row in manifest['files']:
        data = (stage1/row['path']).read_bytes()
        assert len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256'], \
            'stage 1 file differs from its MANIFEST.json: '+row['path']
    work = Path(work)
    shutil.copytree(stage1/'record', work/'record', ignore=shutil.ignore_patterns('__pycache__'))
    for row in manifest['data_files']:
        source = (stage1/row['path']).resolve()
        assert sha(source) == row['sha256'], 'stage 1 data file differs: '+row['path']
        target = work/'record'/row['placed_at']
        assert not target.exists()
        shutil.copyfile(source, target)
    exe = shutil.which(cxx)
    assert exe is not None, 'C++ compiler not found: '+cxx
    tool = work/'toolchain'
    tool.mkdir()
    (tool/'g++').write_text('#!/bin/sh\nexec "'+exe+'" "$@"\n')
    (tool/'g++').chmod(0o755)
    env = {k: v for k, v in os.environ.items() if k != 'PYTHONOPTIMIZE'}
    env.update(PYTHONDONTWRITEBYTECODE='1', OPENBLAS_NUM_THREADS='1',
               PATH=str(tool)+os.pathsep+os.environ.get('PATH', ''))
    return work/'record', env


def run_inner(pkg, work, out, log, cxx='clang++'):
    record, env = prepare(pkg, work, cxx)
    with open(log, 'w') as fh:
        code = subprocess.run([sys.executable, '-B', str(Path(__file__).resolve()), '--inner', '--out', str(out)],
                              cwd=record, env=env, stdout=fh, stderr=subprocess.STDOUT).returncode
    return code


# ------------------------------------------------------------------ inner helpers
def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def inner_function(outer, name, cells):
    """Rebuild a nested function (lambda or def) of `outer` from its code object."""
    codes = [c for c in outer.__code__.co_consts if isinstance(c, types.CodeType) and c.co_name == name]
    assert len(codes) == 1, (outer.__name__, name)
    code = codes[0]
    assert set(code.co_freevars) == set(cells), (code.co_freevars, cells)
    closure = tuple(types.CellType(cells[v]) for v in code.co_freevars)
    return types.FunctionType(code, outer.__globals__, name, None, closure)


class _Identities:
    def __init__(self, degree):
        self.degree = degree

    def verify(self, folder, only=None):
        return dict(status='VERIFIED', polynomial_identity_checks=6, angular_degree=self.degree)


class _Capture:
    """Stub for the lift helper: records the arguments of check()."""
    def __init__(self):
        self.calls = []

    def check(self, core, kind, regions, eta, tents, basis_intervals, budget, probes, **kw):
        self.calls.append(dict(core=core, kind=kind, regions=tuple(regions), eta=eta, tents=tents,
                               basis_intervals=basis_intervals, budget=budget, probes=probes, **kw))
        return dict(status='CAPTURED_NOT_CHECKED')


class _Factor:
    def clear_cache(self):
        pass

    def accept(self, *a, **k):
        raise AssertionError('not used')


def capture_verify(coord, folder):
    snap = {}
    original = coord.load_module

    def wrapped(name, path):
        module = original(name, path)
        if name == 'cancellation_preserving_pair_budget_exact':
            cls = module.CancellationPreservingPairBudget
            init = cls.__init__

            def capturing_init(self, *a, **k):
                init(self, *a, **k)
                snap.update(dict(sys._getframe(1).f_locals))
                snap['refined_budget'] = self
            cls.__init__ = capturing_init
        return module
    coord.load_module = wrapped
    try:
        result = coord.verify(folder, 'pair', 1)
    finally:
        coord.load_module = original
    assert result['status'] == PAIR_STATUS, result['status']
    assert 'refined_budget' in snap
    return result, snap


def main_inner(out):
    started = time.time()
    record = Path.cwd()
    sys.path.insert(0, str(HERE))
    sys.path.insert(0, str((record / COORDINATOR).parent))
    coord = load_module('v7_coordinator_for_kernel_link', record/COORDINATOR)
    folder = (record/CANDIDATE).resolve()
    result, L = capture_verify(coord, folder)
    core = L['core']
    angular = core.gram()[2]
    assert angular == 18 and core.NS == 7 and core.BREAKS[-1] == Q(7, 4)
    capture = _Capture()
    coord._WORKER_STATE = dict(
        parent_ready=True, folder=folder, core=core, budget=L['refined_budget'],
        identities=_Identities(angular), angular_degree=angular,
        expected_jobs=frozenset(L['jobs']), cone=L['cone'], coefficients=L['coefficients'],
        flags=L['flags'], joint_polynomials=L['joint_polynomials'], conditional=L['conditional'],
        tents=L['tents'], basis_intervals=L['basis_intervals'], lift_probes=L['lift_probes'],
        caps=L['caps'], lift_helper=capture, allocations=L['allocations'],
        spherical_rays=L['spherical_rays'], graph_weights=L['graph_weights'],
        outer_graph_weight=L['outer_graph_weight'], narrow_graph_weight=L['narrow_graph_weight'],
        moderate_graph_weight=L['moderate_graph_weight'], factor_backend=_Factor())
    branches = {}
    try:
        for i, j in itertools.combinations_with_replacement(range(NREG), 2):
            job = ('outside', (i, j))
            assert job in L['jobs']
            packet = coord._branch_worker(job)
            call = capture.calls[-1]
            assert call['kind'] == 'outside' and call['regions'] == (i, j)
            assert str(call['eta']) == packet['cyclic']['eta']
            C, dd, domain, report = L['lift_helper'].build_signed(core, 'outside', (i, j), call['eta'],
                                                                  call['spherical_rays'])
            assert report['status'] == 'EXACT_SIGNED_SPHERICAL_BRANCH_IDENTITIES_VERIFIED'
            assert [tuple(map(Q, d)) for d in domain] == [tuple(core.BREAKS[i:i+2]), tuple(core.BREAKS[j:j+2]), (Q(0), Q(1))]
            sparse = {tuple(map(int, idx)): int(v) for idx, v in zip(itertools.product(*map(range, C.shape)), C.flat) if v}
            branches[(i, j)] = dict(C=sparse, dd=int(dd), shape=tuple(C.shape), call=call,
                                    eta=call['eta'], packet=packet['cyclic'])
    finally:
        coord._WORKER_STATE = None
    # every call received the same correction objects
    first = capture.calls[0]
    for call in capture.calls:
        for key in ['tents', 'caps', 'allocations', 'spherical_rays', 'basis_intervals', 'budget']:
            assert call[key] is first[key], key
        for key in ['graph_weights', 'outer_graph_weight', 'narrow_graph_weight', 'moderate_graph_weight']:
            assert call[key] == first[key]
    import stage1_terms
    pos = stage1_terms.build_all(core, first, L, NREG)
    checks = stage1_terms.singleton_checks(core, first, L, NREG, pos)
    payload = dict(
        breaks=[str(b) for b in core.BREAKS], angular_degree=angular,
        branches={k: dict(C=v['C'], dd=v['dd'], shape=v['shape'], eta=str(v['eta']),
                          nested=v['packet']['nested_cyclic_correction'],
                          joint=v['packet']['conditional_polynomial_cyclic_correction'],
                          mixed=v['packet']['conditional_cyclic_correction'])
                  for k, v in branches.items()},
        corrections={cell: sym.key() for cell, sym in pos['cells'].items()},
        endpoints={cell: sym.key() for cell, sym in pos['endpoints'].items()},
        components={cell: {name: s.key() for name, s in comp.items()} for cell, comp in pos['components'].items()},
        singleton_checks=checks, coordinator_status=result['status'],
        pair_pieces=len(result['pairs']), seconds=time.time()-started)
    Path(out).write_bytes(pickle.dumps(payload))
    print('EXTRACT_STAGE1_DONE', len(branches), 'branches', len(pos['cells']), 'cells',
          checks['points'], 'singleton points', '%.1fs' % (time.time()-started), flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--inner', action='store_true')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    assert a.inner
    main_inner(a.out)
