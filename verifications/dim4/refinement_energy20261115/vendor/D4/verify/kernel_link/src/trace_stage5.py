"""Trace stage 5's kernel code exactly on every region-pair cell of B_5.

Loads the unmodified src/kcert.py and src/kformula.py of stage 5 from --stage5,
replaces their number type by symalg.Sym (exact rational functions in r, s, t
plus kink atoms) and runs

  kcert.minusK_lb(r, s, t, ra, rb)      the function whose lower bounds stage 5 certifies
  kformula.pieces(r, s, t, ra, rb)      stage 5's explicit formula (smooth + pos)

with symbolic r, s, t, and with ra, rb the decision intervals kcert uses on the
root boxes of the cell. Patched names and why:

kcert:  T -> Sym; A, arb -> exact Fraction; ZERO -> 0. The kink lower-bound helpers
        relu_lb, clip01_lb, ind_lb, capacity_lb, graph_lb are replaced by the EXACT
        kink (symalg.relu, clip01, ind_le0, capacity, band). kcert calls each helper on
        an argument computed by its own code; the trace records that argument. Each
        replaced helper returns, on a box, a constant lower bound of the exact kink or
        the exact branch. The written lower-envelope proof is in kernel_box_method.md
        section 3; this trace checks their exact arguments, not the interval primitives.
kformula: A, arb -> exact; relu, clip01, capacity, cap_prob, cap_prob_self, graph_pos ->
        exact kinks (these helpers return Arb enclosures [0,1] when undecided, which
        carry no symbolic information). Everything else runs as written.

Output: a pickle with, for each ordered cell (i, j), the canonical keys of -K
from kcert and from kformula, and the log of kcert's D.lower() calls.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import importlib
import pickle
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import symalg as sa

AA = Q(133, 250)
DELTA = Q(1, 2 ** 40)


def graph_exact(kf, ra, rb, t):
    """kcert.graph_lb with the exact band indicator in place of its lower bound."""
    total = sa.Sym.const(0)
    for w, lr, hr, lt, ht in kf.GRAPH_SPECS:
        inr = lambda R: (lr <= R[0] and R[1] <= hr, R[1] < lr or R[0] > hr)
        ya, na = inr(ra)
        yb, nb = inr(rb)
        if na or nb:
            continue
        assert ya and yb, 'box straddles a graph radius'
        total = total + sa.band(t, lt, ht) * w
    return total


def cap_prob_exact(row, r, t):
    """The exact cap probability as written in kformula.cap_prob (exact kinks)."""
    d, a2, sg = Q(row['d']), Q(row['a2']), Q(row['sigma'])
    if sg * sg == 1:
        return sa.ind_le0(r * r + (d * d - a2) - 2 * (sg * d) * r * t)
    c = (r * r + (d * d - a2)) / (2 * d * r)
    if t == -1:
        return sa.ind_le0(c + sg)
    D = (1 - sg * sg) * (1 - t * t)
    return sa.clip01((1 - (c - sg * t) / D.sqrt()) / 2)


def cap_prob_self_exact(row, r):
    d, a2, sg = Q(row['d']), Q(row['a2']), Q(row['sigma'])
    return sa.ind_le0(r * r - 2 * (sg * d) * r + (d * d - a2))


def load(stage5_src):
    stage5_src = Path(stage5_src).resolve()
    sys.path.insert(0, str(stage5_src))
    kf = importlib.import_module('kformula')
    kc = importlib.import_module('kcert')
    assert Path(kf.__file__).resolve().parent == stage5_src
    assert Path(kc.__file__).resolve().parent == stage5_src
    assert kc.kf is kf
    # kcert: exact number type, exact kinks.
    kc.T = sa.Sym
    kc.A = lambda q: Q(q)
    kc.arb = lambda x: Q(x)
    kc.ZERO = Q(0)
    kc.relu_lb = sa.relu
    kc.clip01_lb = sa.clip01
    kc.ind_lb = sa.ind_le0
    kc.capacity_lb = sa.capacity
    kc.graph_lb = lambda ra, rb, t: graph_exact(kf, ra, rb, t)
    kc._RW.clear()
    kc._GA.clear()
    return kf, kc


def patch_kformula(kf):
    kf.A = lambda q: Q(q)
    kf.arb = lambda x: Q(x)
    kf.relu = sa.relu
    kf.clip01 = sa.clip01
    kf.capacity = sa.capacity
    kf.cap_prob = cap_prob_exact
    kf.cap_prob_self = cap_prob_self_exact
    kf.graph_pos = lambda ra, rb, t: graph_exact(kf, ra, rb, t)


def cells(kf, j=5):
    R = kf.BREAKS
    return [(i, k) for i in range(j) for k in range(j)]


def trace(stage5_src, out):
    t0 = time.time()
    kf, kc = load(stage5_src)
    R = kf.BREAKS
    assert len(R) == 8 and R[0] == 1
    r, s, t = sa.var('r'), sa.var('s'), sa.var('t')
    result = dict(breaks=[str(x) for x in R], a=str(AA), cells={}, lower_log={}, endpoints={})
    # kcert first (kformula patches below change shared kf names only after this loop).
    for i, k in cells(kf):
        box = ((R[i], R[i + 1]), (R[k], R[k + 1]), (Q(-1), -AA))
        sa.Sym.BOX = box
        sa.Sym.GENERIC_T = (Q(-1) + DELTA, -AA)
        sa.Sym.LOWER_LOG = []
        ra, rb = kc.dec(box[0], i > 0), kc.dec(box[1], k > 0)
        val = kc.minusK_lb(r, s, t, ra, rb)
        assert isinstance(val, sa.Sym)
        result['cells'][(i, k)] = dict(kcert=val.key(), ra=[str(x) for x in ra], rb=[str(x) for x in rb])
        result['lower_log'][(i, k)] = list(sa.Sym.LOWER_LOG)
    patch_kformula(kf)
    for i, k in cells(kf):
        box = ((R[i], R[i + 1]), (R[k], R[k + 1]), (Q(-1), -AA))
        sa.Sym.BOX = box
        ra, rb = kc.dec(box[0], i > 0), kc.dec(box[1], k > 0)
        smooth, pos = kf.pieces(r, s, t, ra, rb)
        val = sa.Sym._wrap(smooth) + sa.Sym._wrap(pos)
        result['cells'][(i, k)]['kformula'] = val.key()
    for i, k in cells(kf):
        box = ((R[i], R[i + 1]), (R[k], R[k + 1]), (Q(-1), Q(-1)))
        sa.Sym.BOX = box
        ra, rb = kc.dec(box[0], i > 0), kc.dec(box[1], k > 0)
        t = sa.Sym.const(-1)
        smooth, pos = kf.pieces(r, s, t, ra, rb)
        exact = sa.Sym._wrap(smooth) + sa.Sym._wrap(pos)
        lower = kc.minusK_lb(r, s, t, ra, rb)
        # At t=-1 the Taylor code drops the nonnegative caps with |sigma|<1.
        # Reconstruct precisely those dropped terms from the same data. All
        # coefficients are nonnegative Bernstein coefficients on these cells.
        dropped = sa.Sym.const(0)
        for row in kf.origin['pointed_caps']:
            if abs(Q(row['sigma'])) == 1:
                continue
            assert Q(row['weight']) >= 0 and all(Q(w) >= 0 for w in row['w'])
            wa = kf.radial_weight(row, r, i, ra[1])
            wb = kf.radial_weight(row, s, k, rb[1])
            dropped += Q(row['weight']) * (wb*cap_prob_exact(row,r,t)+wa*cap_prob_exact(row,s,t))/2
        assert (exact-lower-dropped).is_zero(), 'endpoint omitted-cap identity'
        result['endpoints'][i,k] = dict(exact=exact.key(), lower=lower.key(), dropped=dropped.key())
    result['seconds'] = time.time() - t0
    Path(out).write_bytes(pickle.dumps(result))
    return result


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage5', required=True, help='stage5_kernelbox/src folder')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    res = trace(a.stage5, a.out)
    print('TRACE_STAGE5_DONE', len(res['cells']), 'cells', '%.1fs' % res['seconds'])
