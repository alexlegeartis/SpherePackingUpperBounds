"""Coordinator spot checks for finding 1114 (independent of kformula/kcert):
(1) additivity: raw_v7 E on 3- and 4-point stars equals A + 2 sum h + sum K via polarization;
(2) -K - lam*q >= 0 at rational points of B_1 near the dip and at the band corners, using kernel_poly enclosures.
Run from any folder. Writes paths.OUT/coord_check_kernelbox.json."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util, sys, random, json
sys.set_int_max_str_digits(0)
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import paths
spec = importlib.util.spec_from_file_location("kernel_poly", HERE / "kernel_poly.py")
kp = importlib.util.module_from_spec(spec); spec.loader.exec_module(kp)
raw = kp.raw; A = kp.A
O = (Q(0),)*4
def dot(x,y): return sum(a*b for a,b in zip(x,y))
def sep_ok(pts):
    allp = [O]+pts
    return all(dot(tuple(a-b for a,b in zip(p,q)), tuple(a-b for a,b in zip(p,q))) >= 1 for i,p in enumerate(allp) for q in allp[i+1:])
random.seed(7)
# (1) additivity on rational stars of 3 and 4 points inside radius 7/4
out = {"additivity": [], "band_points": []}
def rand_pt():
    while True:
        v = tuple(Q(random.randint(-14,14),8) for _ in range(4))
        n = dot(v,v)
        if 1 <= n <= Q(49,16): return v
tests = 0
while tests < 4:
    k = 3 if tests < 2 else 4
    pts = [rand_pt() for _ in range(k)]
    if not sep_ok(pts): continue
    E = raw.evaluate([O]+pts, O)["total"]
    h = [raw.evaluate([O,p], O)["total"] for p in pts]  # singleton energies: A + 2h(p) + K(p,p)
    # pairwise reconstruction: E = A + sum_p (E1(p)-A) + sum_{p!=q} K(p,q)
    lo = A + sum(x[0]-A for x in h); hi = A + sum(x[1]-A for x in h)
    for i,p in enumerate(pts):
        for q in pts[i+1:]:
            kk = kp.K(p,q); lo += 2*kk[0]; hi += 2*kk[1]
    ok = not (hi < E[0] or lo > E[1])
    out["additivity"].append({"k": k, "E": [float(E[0]), float(E[1])], "recon": [float(lo), float(hi)], "overlap": ok})
    tests += 1
# (2) band points: radii r,s in [1,R1], t in [-1,-0.532]
R1 = Q(501473,500000); a = Q(133,250); lam1 = Q(3,10**6)
q = lambda t: t*t*(t*t-a*a)
def vec(r, t, s):
    # y = r e1, z = s (t e1 + sqrt(1-t^2) e2) needs irrational; use rational t with rational sqrt: pick t = -(p^2-q^2)/(p^2+q^2)
    return None
pts = []
# rational cos values with rational sine: t = -(m^2-n^2)/(m^2+n^2)
cands = []
for m in range(1,30):
    for n in range(0,m):
        t = -Q(m*m-n*n, m*m+n*n)
        if -1 <= t <= -a: cands.append((t, Q(2*m*n, m*m+n*n)))
cands = sorted(set(cands))
sel = [c for c in cands if -0.9 < c[0] < -0.83][:12] + [cands[0], cands[-1]] + random.sample(cands, 8)
for (t, sn) in sel:
    for (r, s) in [(Q(1),Q(1)), (Q(1),R1), (R1,R1), (Q(1000500,1000000), Q(1))]:
        y = (r,0,0,0); z = (s*t, s*sn, 0, 0)
        if dot(tuple(a1-b1 for a1,b1 in zip(y,z)), tuple(a1-b1 for a1,b1 in zip(y,z))) < 1: continue
        kk = kp.K(y,z)
        val_lo = -kk[1] - lam1*q(t)
        out["band_points"].append({"t": float(t), "r": float(r), "s": float(s), "negK_lo": float(-kk[1]), "margin_lo": float(val_lo), "ok": val_lo >= 0})
bad = [b for b in out["band_points"] if not b["ok"]]
out["summary"] = {"additivity_all_overlap": all(x["overlap"] for x in out["additivity"]), "band_points": len(out["band_points"]), "band_failures": len(bad)}
print(json.dumps(out["summary"]))
print(json.dumps(out["additivity"], indent=0))
print(json.dumps(sorted(out["band_points"], key=lambda b: b["margin_lo"])[:6], indent=0))
paths.OUT.mkdir(parents=True, exist_ok=True)
(paths.OUT / "coord_check_kernelbox.json").write_text(json.dumps(out, indent=1))
