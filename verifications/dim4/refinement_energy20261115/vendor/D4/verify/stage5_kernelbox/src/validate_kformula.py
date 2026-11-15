"""Check kformula against the retained evaluator (kernel_poly.K, raw_v7.unlifted_pair).

Points: y = (r,0,0,0), z = rho (c, s, 0, 0), c = -(1-u^2)/(1+u^2), s = 2u/(1+u^2),
so t = c is rational. Pass condition: the kformula ball meets the evaluator's
rational enclosure for every point (both are enclosures of the same exact value).
Run from any folder. Writes paths.OUT/validate_kformula.json.
"""
from fractions import Fraction as Q
from pathlib import Path
import json, sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths
import kformula as kf
import kernel_poly as kp
from flint import arb

R = kf.BREAKS
radii = [Q(1), (Q(1) + R[1]) / 2, R[1], (R[1] + R[2]) / 2, R[2], R[3], (R[3]+R[4])/2, R[4], R[5], Q(103,100)]
us = [Q(k, 40) for k in (1, 4, 8, 11, 13, 16, 20, 21)]
pairs = [(radii[0], radii[0]), (radii[1], radii[2]), (radii[2], radii[2]), (radii[0], radii[3]),
         (radii[4], radii[5]), (radii[6], radii[7]), (radii[8], radii[8]), (radii[0], radii[8]), (radii[9], radii[2])]
rows, worst = [], 0
for r, rho in pairs:
    for u in us:
        c = -(1 - u * u) / (1 + u * u); s = 2 * u / (1 + u * u)
        y = (r, Q(0), Q(0), Q(0)); z = (rho * c, rho * s, Q(0), Q(0))
        lo, hi = kp.K(y, z)
        mine = kf.minusK_point(r, rho, c)
        ref = kf.A(-hi).union(kf.A(-lo))
        ok = mine.overlaps(ref)
        rows.append(dict(r=str(r), rho=str(rho), t=str(c), ref_lo=float(-hi), ref_hi=float(-lo), kformula=float(mine.mid()), ok=bool(ok)))
        worst = max(worst, abs(float(mine.mid()) + float((lo + hi) / 2)))
        assert ok, rows[-1]
q0rows = []
for r in radii[:9]:
    lo, hi = kf_ref = kp.raw.unlifted_pair(r * r)
    mine = kf.q0(kf.A(r), (r, r))
    q0rows.append(dict(r=str(r), ref=[float(lo), float(hi)], kformula=float(mine.mid())))
    # unlifted_pair clips the upper endpoint at 0; the true value lies in [lo, hi]
    assert mine.overlaps(kf.A(lo).union(kf.A(hi))), q0rows[-1]
out = dict(status="KFORMULA_MATCHES_RAW_V7", points=len(rows), max_mid_difference=worst, rows=rows, q0=q0rows)
paths.OUT.mkdir(parents=True, exist_ok=True)
(paths.OUT / "validate_kformula.json").write_text(json.dumps(out, indent=1))
print(out["status"], len(rows), worst)
