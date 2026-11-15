"""Bonus: compare kernel.py with stage 5's arb kernel formula (kformula.py)
on the kernel-box domain [1, R5]^2 x [-1, -133/250], at sample points.

This is a sampled comparison, not a proof: at each rational point (r, s, t)
it checks that kernel.py's value of -K(x, y) (exact rational, or an arb ball
for the square-root and exp parts) and stage 5's arb ball for -K overlap.
Stage 5's balls have radius up to about 1e-4 (its q0 enclosure), so this
test is weak. A sharper test uses stage 5's own Gram, spherical, count and
conditional pieces (tight balls, radius < 1e-60): kernel.py's exact
polynomial part must lie inside.

    python3 bonus_stage5.py --out FILE [--data DIR]
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import itertools
import json
import sys

HERE = Path(__file__).resolve().parent
PKG = HERE.parents[1]


def overlap(a, b):
    d = a-b
    return d.contains(0) if hasattr(d, 'contains') else bool(0 in d)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--data', type=Path, default=PKG/'data'/'v7')
    args = parser.parse_args()
    sys.path.insert(0, str(HERE))
    import kernel as kmod
    import flint
    K = kmod.Kernel(args.data)
    sys.path.insert(0, str(PKG/'verify'/'stage5_kernelbox'/'src'))
    import kformula
    flint.ctx.prec = kmod.PREC
    R5 = K.breaks[5]
    # radii in (B_p, B_p+1] for p = 0..4, plus r = 1; cosines in [-1, -133/250]
    radii = [Q(1)]
    for p in range(5):
        lo, hi = K.breaks[p:p+2]
        radii += [lo+(hi-lo)/3, hi]
    cosines = [Q(-1), Q(-9, 10), Q(-3, 4), Q(-7, 10), Q(-3, 5), Q(-133, 250)]
    rows, ok, n, widest, poly_tight, finite = [], True, 0, 0.0, 0, 0
    for r, s, t in itertools.product(radii, radii, cosines):
        if r > s:
            continue
        assert r <= R5 and s <= R5
        mine = -K.pair_kernel_point(r, s, t)
        mine = kmod._add(mine, K.caps_k(r, s, t))
        mine = kmod._add(mine, K.tents_k(r, s, t))
        mine = kmod._add(mine, K.allocs_k(r, s, t))
        mine = kmod._add(mine, K.graphs_k(r, s, t))
        mine = kmod._add(mine, K.gauge_k(r, s))
        theirs = kformula.minusK_point(r, s, t)
        good = overlap(kmod.to_arb(mine), theirs)
        # Sharper: the polynomial components alone, built from stage 5's own
        # pieces (Gram, spherical, count and conditional terms), whose balls
        # are tight. kernel.py's exact value must lie inside.
        A = kformula.A
        ra, rb = (r, r), (s, s)
        fa = kformula.features(A(r), kformula.region_of(*ra))
        fb = kformula.features(A(s), kformula.region_of(*rb))
        Z = kformula.zonal_list(A(t), len(kformula.GRAM)-1)
        poly5 = sum((Z[d]*kformula.quad_form(kformula.GRAM[d], fa, fb)
                     for d in range(len(kformula.GRAM))), flint.arb(0))
        for row in kformula.origin['appended_spherical_local_energy_rays']:
            R = Q(row['R'])
            if kformula.ind_le(*ra, R) and kformula.ind_le(*rb, R):
                poly5 -= A(row['weight'])*sum((A(c)*Z[d] for d, c in enumerate(row['scalar_coefficients'])),
                                              flint.arb(0))
        poly5 += A(kformula.K_count(ra, rb)+kformula.K_cond(ra, rb))
        exact = K.pair_kernel_point(r, s, t)
        tight = overlap(kmod.to_arb(exact), poly5) and float(poly5.rad()) < 1e-60
        poly_tight += int(tight)
        good = good and tight
        n += 1
        if theirs.is_finite():
            finite += 1
            widest = max(widest, float(theirs.rad()))
        if not good:
            ok = False
            rows.append(dict(r=str(r), s=str(s), t=str(t), mine=str(mine), stage5=str(theirs)))
    result = dict(status='PASS' if ok else 'FAIL', points=n, failures=rows[:10],
                  points_with_finite_stage5_ball=finite,
                  widest_finite_stage5_ball_radius=widest,
                  polynomial_part_points_inside_tight_stage5_ball=poly_tight,
                  scope='sampled overlap of arb balls; not a proof')
    args.out.write_text(json.dumps(result, indent=1)+'\n')
    print(json.dumps(result)[:1000])
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
