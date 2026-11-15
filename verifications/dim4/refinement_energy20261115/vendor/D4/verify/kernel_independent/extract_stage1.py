"""Run the stage 1 V7 coordinator in a scratch copy and compare its exact
objects with the independent kernel of kernel.py.

Usage (normally called by run.py, with cwd = the scratch record copy):

    python3 extract_stage1.py --record REC --data DATA --out OUT [--workers N]

REC is a scratch copy of verify/stage1_v7/record with the four data/v7 files
placed as stage 1's run.py places them. DATA is the data/v7 folder read by
kernel.py. OUT receives one JSON file per compared object.

How stage 1's objects are obtained. The coordinator's verify(stage='cyclic')
runs unchanged through every parent gate: graph capacities, count and
spherical replays, origin identity, objective and the 7 fresh pair proofs.
It then forks workers that call _branch_worker(job) for each of the 112
cases. This script replaces only the worker function. The replacement calls
the original worker with two substitutions inside the forked child:
  - the per-case identity audit is skipped (it does not change any object);
  - lift_helper.check, the cyclic subdivision prover, is replaced by a
    capture that calls lift_helper.build_signed(core, kind, regions, eta,
    spherical_rays) with the exact arguments the original check receives.
    That is the first statement of the original check, and it returns the
    numerator C and denominator dd whose sign the original check proves.
The capture then compares C/dd with kernel.py's polynomial and evaluates the
nonpolynomial corrections at sample points. The worker returns an UNKNOWN
packet, so the coordinator never reports a certificate. Its final coverage
assertion then fails by design; this script checks its own 112 outputs.
No stage 1 file is modified.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import importlib.util
import itertools
import json
from numbers import Integral
import os
import sys
import time
import traceback

CANDIDATE = 'current/candidates/round9_six_graphs_degree18'
COORDINATOR = 'current/kernel_integration_route/verify_combined_binary_six_exterior_graphs_small_cap_bounds.py'


class Captured(Exception):
    pass


def atomic_json(path, value):
    tmp = Path(str(path)+'.tmp')
    tmp.write_text(json.dumps(value, indent=1)+'\n')
    tmp.replace(path)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def canonical_hash(terms):
    h = hashlib.sha256()
    for key in sorted(terms):
        value = terms[key]
        h.update(('%s:%d/%d\n' % (','.join(map(str, key)), value.numerator, value.denominator)).encode())
    return h.hexdigest()


# ----------------------------------------------------------------------
# comparison helpers (kernel side only uses kernel.py)

def compare_polynomial(C, dd, mine):
    """Exact coefficientwise comparison of stage 1 C/dd with an fmpq_mpoly."""
    import numpy as np
    if isinstance(dd, bool) or not isinstance(dd, Integral) or dd <= 0:
        raise ValueError('stage 1 polynomial denominator must be a positive integer')
    if C.ndim != 3 or any(n <= 0 for n in C.shape):
        raise ValueError('stage 1 polynomial tensor must have three nonempty axes')
    stage = {}
    for idx in zip(*np.nonzero(C)):
        idx = tuple(int(i) for i in idx)
        coefficient = C[idx]
        if isinstance(coefficient, bool) or not isinstance(coefficient, Integral):
            raise ValueError('stage 1 polynomial coefficient must be an integer')
        stage[idx] = Q(int(coefficient), int(dd))
    own = {}
    for key, value in mine.to_dict().items():
        own[tuple(int(k) for k in key)] = Q(int(value.p), int(value.q))
    mismatches = []
    for key in sorted(set(stage) | set(own)):
        a, b = stage.get(key, Q(0)), own.get(key, Q(0))
        if a != b:
            mismatches.append(dict(monomial=list(key), stage1=str(a), independent=str(b)))
            if len(mismatches) >= 5:
                break
    degrees = [max((k[i] for k in stage), default=0) for i in range(3)]
    return dict(equal=not mismatches, stage1_terms=len(stage), independent_terms=len(own),
                max_degrees=degrees, shape=list(C.shape), mismatches=mismatches,
                stage1_sha256=canonical_hash(stage) if not mismatches else None)


def within(enclosure, mine, flint):
    """Is my value inside stage 1's closed [lo, hi]? Exact if mine is rational."""
    lo, hi = map(Q, enclosure)
    if not lo <= hi:
        return False, 'bad enclosure'
    if isinstance(mine, Q):
        ok = lo <= mine <= hi
        return ok, ('exact_equal' if ok and lo == mine == hi else ('exact_inside' if ok else 'exact_outside'))
    a = mine-flint.arb(flint.fmpq(lo.numerator, lo.denominator))
    b = flint.arb(flint.fmpq(hi.numerator, hi.denominator))-mine
    if a >= 0 and b >= 0:
        return True, 'ball_inside'
    if a < 0 or b < 0:
        return False, 'ball_outside'
    return False, 'undecided'


def sample_coordinates(lo, hi, first_piece):
    """Rational sample positions in (lo, hi], plus lo on the first piece."""
    pts = [hi, lo+(hi-lo)/2, lo+(hi-lo)/7]
    if first_piece:
        pts.append(lo)
    return pts


def sample_points(K, kind, regions):
    pieces = [tuple(K.breaks[p:p+2]) for p in regions]
    axes = [sample_coordinates(lo, hi, p == 0) for (lo, hi), p in zip(pieces, regions)]
    if kind == 'outside':
        axes.append([Q(0), Q(1, 3), Q(5, 7), Q(1)])
    return [tuple(p) for p in itertools.product(*axes)]


def correction_triplets(K, kind, point):
    """The (r, s, t) arguments of the K(x,y) terms of Gamma at a point."""
    if kind == 'full':
        a, b, c = point
        cos = lambda u, v, w: (u*u+v*v-w*w)/(2*u*v)
        return [(a, b, cos(a, b, c)), (a, c, cos(a, c, b)), (b, c, cos(b, c, a))]
    a, b, z = point
    return [(a, b, K.outside_cosine(a, b, z))]


# ----------------------------------------------------------------------

_CTX = {}


def worker(job):
    """Replacement for the coordinator's _branch_worker (module level, so it pickles)."""
    K, coord, out, only, flint = (_CTX[k] for k in ('K', 'coord', 'out', 'only', 'flint'))
    original_worker = _CTX['original_worker']
    state = coord._WORKER_STATE
    kind, regions = job
    name = '%s_%s.json' % (kind, '_'.join(map(str, regions)))
    started = time.time()
    report = dict(kind=kind, regions=list(regions), status='ERROR')
    try:
        if only is not None and job not in only:
            report['status'] = 'SKIPPED'
            raise Captured
        if job == ('full', (0, 0, 0)):
            import contracts
            atomic_json(out/'contracts.json', contracts.compare(K, state))
            atomic_json(out/'pair.json', compare_pair(K, state, flint))
        lift = state['lift_helper']
        captured = {}

        def capture(core, kind_, regions_, eta, tents, basis_intervals, budget, probes, **kw):
            assert (kind_, tuple(regions_)) == job
            C, dd, domain, identity = lift.build_signed(core, kind_, regions_, eta, kw['spherical_rays'])
            captured.update(C=C, dd=dd, domain=domain, eta=eta, identity=identity, core=core,
                            tents=tents, basis_intervals=basis_intervals, budget=budget,
                            probes=probes, kw=kw)
            raise Captured

        class SkipIdentities:
            @staticmethod
            def verify(folder_, only=None):
                return dict(status='VERIFIED', polynomial_identity_checks=6,
                            angular_degree=state['angular_degree'])

        lift.check = capture
        state['identities'] = SkipIdentities
        try:
            original_worker(job)
        except Captured:
            pass
        assert captured, 'the original worker did not reach the cyclic check'
        mine = K.gamma_poly(kind, regions)
        poly = compare_polynomial(captured['C'], captured['dd'], mine)
        del captured['C']
        mydomain = [tuple(K.breaks[p:p+2]) for p in regions]
        if kind == 'outside':
            mydomain.append((Q(0), Q(1)))
        domain_ok = [tuple(map(Q, d)) for d in captured['domain']] == mydomain
        samples = sample_corrections(K, kind, regions, captured, lift, flint)
        report.update(status='PASS' if poly['equal'] and domain_ok and samples['ok'] else 'FAIL',
                      polynomial=poly, domain_equal=domain_ok, eta=str(captured['eta']),
                      stage1_signed_spherical_identity=captured['identity']['status'],
                      samples=samples)
    except Captured:
        pass
    except BaseException as error:
        report.update(status='ERROR', error='%s: %s' % (type(error).__name__, error),
                      traceback=traceback.format_exc())
    report['seconds'] = round(time.time()-started, 2)
    report['pid'] = os.getpid()
    if report['status'] != 'SKIPPED':
        atomic_json(out/name, report)
    print('KERNEL_INDEPENDENT', kind, regions, report['status'], report['seconds'], flush=True)
    return dict(kind=kind, regions=list(regions), identities={},
                cyclic=dict(status='UNKNOWN', kind=kind, regions=list(regions),
                            note='kernel_independent extraction; not a proof'))



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--record', type=Path, required=True)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--kernel-dir', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--only', default=None, help='comma list like full:0,0,0;outside:5,5 (testing)')
    parser.add_argument('--perturb', default=None, help='negative control in kernel.py')
    args = parser.parse_args()
    if not __debug__:
        raise SystemExit('run without -O')
    record = args.record.resolve()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(args.kernel_dir.resolve()))
    import kernel as kmod
    import flint
    perturb = None
    if args.perturb == 'implementation':
        def perturb(K):
            # Negative control: change one spherical ray weight by 1e-9 relative.
            K.sph[0]['weight'] *= 1+Q(1, 10**9)
    K = kmod.Kernel(args.data.resolve(), perturb=perturb)
    only = None
    if args.only:
        only = set()
        for item in args.only.split(';'):
            kind, regs = item.split(':')
            only.add((kind, tuple(int(x) for x in regs.split(','))))

    # As when the coordinator runs as a script: its folder is sys.path[0].
    sys.path.insert(0, str((record/COORDINATOR).parent))
    coord = load('stage1_v7_coordinator', record/COORDINATOR)
    original_worker = coord._branch_worker
    folder = record/CANDIDATE

    _CTX.update(K=K, coord=coord, out=out, only=only, flint=flint,
                original_worker=original_worker)
    coord._branch_worker = worker
    started = time.time()
    error = None
    try:
        coord.verify(folder, stage='cyclic', workers=args.workers)
    except BaseException as exc:
        error = '%s: %s' % (type(exc).__name__, exc)
    # The coordinator's progress file carries its exact objective (U_raw).
    progress = sorted(folder.glob('*CYCLIC_PROBE.json'))
    summary = dict(coordinator_exception=error, seconds=round(time.time()-started, 1))
    if progress:
        data = json.loads(progress[-1].read_text())
        summary.update(stage1_objective=data.get('objective'),
                       stage1_base_objective=data.get('base_objective'),
                       stage1_origin_identity=data.get('origin_identity'),
                       stage1_pairs_certified=sum(p.get('status') == 'CERTIFIED' for p in data.get('pairs', [])),
                       stage1_unproved_jobs=len(data.get('parallel', {}).get('unproved_jobs', [])))
    atomic_json(out/'coordinator.json', summary)
    print('KERNEL_INDEPENDENT_EXTRACTION_DONE', json.dumps(summary)[:500], flush=True)


def compare_pair(K, state, flint):
    """Pair term on the 7 pieces: exact polynomial parts, then the cap self term."""
    budget = state['budget']
    base = budget.base
    core = state['core']
    caps = state['caps']
    basis_intervals = state['basis_intervals']
    pair_module = sys.modules[type(base).__module__]
    diagonal_bounds = pair_module.diagonal_bounds
    rows = []
    ok = True
    Fpoly = K.F_of_alpha_r2()
    for p in range(K.np):
        fp, dp = base.polynomials[p]
        mine = K.pair_poly(p)
        width = max(len(dp), len(mine))
        pad = lambda v: list(v)+[Q(0)]*(width-len(list(v)))
        d_equal = pad([Q(x) for x in dp]) == pad(mine)
        f_equal = list(map(Q, fp)) == Fpoly+[Q(0)]*(len(fp)-len(Fpoly))
        # Cap self term J: piecewise polynomial of degree <= 6 on cells whose
        # ends are the pivot cutoffs and the roots of r^2 - 2 d sigma r + d^2 - a^2.
        lo, hi = K.breaks[p:p+2]
        cells = J_cells(K, lo, hi)
        points = []
        for a, b in cells:
            points.extend(a+(b-a)*Q(k, 9) for k in range(1, 9))
        points.extend(e for e in [lo, hi] if p == 0 or e == hi)
        points.extend(cap['cutoff'] for cap in K.caps
                      if cap['cutoff'] is not None and lo < cap['cutoff'] < hi)
        bad = []
        for r in points:
            s1 = diagonal_bounds((r, r), p, caps, core, basis_intervals)
            mine_J = K.J_equivalent_point(r)
            if not (Q(s1[0]) == Q(s1[1]) == mine_J):
                bad.append(dict(r=str(r), stage1=list(map(str, s1)), independent=str(mine_J)))
        row = dict(piece=p, interval=[str(lo), str(hi)], D_equal=d_equal, F_equal=f_equal,
                   D_degree=len(mine)-1, J_cells=len(cells), J_points=len(points),
                   J_mismatches=bad[:5])
        ok = ok and d_equal and f_equal and not bad
        rows.append(row)
    return dict(status='PASS' if ok else 'FAIL', pieces=rows)


def J_cells(K, lo, hi):
    """Open cells of (lo, hi) on which every cap self indicator and cutoff is constant.
    Irrational break points are bracketed by rationals 2^-80 apart; cells are
    shrunk to rational endpoints strictly inside."""
    import math
    cuts = []
    for cap in K.caps:
        if cap['cutoff'] is not None and lo < cap['cutoff'] < hi:
            cuts.append((cap['cutoff'], cap['cutoff']))
        d, a2, sigma = cap['d'], cap['a2'], cap['sigma']
        disc = d*d*sigma*sigma-d*d+a2  # roots of r^2 - 2 d sigma r + d^2 - a2
        if disc < 0:
            continue
        scale = 1 << 80
        root = math.isqrt(disc.numerator*scale*scale//disc.denominator)
        for sign in (1, -1):
            r1 = d*sigma+sign*Q(root, scale)
            r2 = d*sigma+sign*Q(root+1, scale)
            a, b = min(r1, r2), max(r1, r2)
            if b > lo and a < hi:
                cuts.append((a, b))
    cuts.sort()
    cells, left = [], lo
    for a, b in cuts:
        if a > left:
            cells.append((left, a))
        left = max(left, b)
    if left < hi:
        cells.append((left, hi))
    return cells


def sample_corrections(K, kind, regions, captured, lift, flint):
    """Compare each nonpolynomial correction at sample points (not a proof)."""
    core = captured['core']
    kw = captured['kw']
    caps = kw['caps']
    allocations = kw['allocations']
    graph_weights = kw['graph_weights']
    tents = captured['tents']
    basis_intervals = captured['basis_intervals']
    budget = captured['budget']
    probes = captured['probes']
    stats = {}
    failures = []

    def record(name, enclosure, mine, point):
        ok, how = within(enclosure, mine, flint)
        stats.setdefault(name, {}).setdefault(how, 0)
        stats[name][how] += 1
        if not ok and len(failures) < 10:
            failures.append(dict(component=name, point=list(map(str, point)), how=how,
                                 stage1=list(map(str, enclosure)), independent=str(mine)))

    for point in sample_points(K, kind, regions):
        single_bounds = [(x, x) for x in point]
        trip = correction_triplets(K, kind, point)
        add = kmod_add
        mine = dict(caps=Q(0), tents=Q(0), allocations=Q(0), gauge=Q(0), graphs=Q(0))
        for r, s, t in trip:
            mine['caps'] = add(mine['caps'], K.caps_k(r, s, t))
            mine['tents'] = add(mine['tents'], K.tents_k(r, s, t))
            mine['allocations'] = add(mine['allocations'], K.allocs_k(r, s, t))
            mine['gauge'] = add(mine['gauge'], K.gauge_k(r, s))
            mine['graphs'] = add(mine['graphs'], K.graphs_k(r, s, t))
        record('caps', lift.cap_bounds(kind, single_bounds, regions, caps, core, basis_intervals),
               mine['caps'], point)
        record('tents', lift.tent_bounds(kind, single_bounds, regions, tents, core, basis_intervals),
               mine['tents'], point)
        record('allocations', lift.allocation_bounds(kind, single_bounds, regions, allocations, core,
                                                      basis_intervals), mine['allocations'], point)
        record('gauge', lift.lift_bounds(kind, single_bounds, regions, probes, budget),
               mine['gauge'], point)
        g = (lift.cyclic_decrement_lower(kind, single_bounds, graph_weights)
             + lift.outer_graph_lower(kind, single_bounds, kw['outer_graph_weight'])
             + lift.narrow_graph_lower(kind, single_bounds, kw['narrow_graph_weight'])
             + lift.moderate_graph_lower(kind, single_bounds, kw['moderate_graph_weight']))
        record('graphs', (g, g), mine['graphs'], point)
    return dict(ok=not failures, points=len(sample_points(K, kind, regions)), outcomes=stats,
                failures=failures)


def kmod_add(a, b):
    import kernel
    return kernel._add(a, b)


if __name__ == '__main__':
    main()
