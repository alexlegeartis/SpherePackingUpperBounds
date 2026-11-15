"""Lemma 1044 (kernelbox) final checker. Fail-closed; exit 0 only if every gate passes.

Gates
 1 exact identities: q = U4/16 + (3/16 - a^2/4) U2 + (1/8 - a^2/4) U0 with nonnegative
   coefficients; s_j = 1 - 1/(2 R_j^2) <= a for j = 1..5; 2 + 2a > 49/16 (far pairs).
 2 validate_kformula.py (kformula vs raw_v7 at 72 points, q0 at 9 radii) exits 0.
 3 kcert point boxes reproduce kformula at 6 points (lower bound <= value, within 1e-30).
 4 replay of the saved split trees kbox_j1.json (lam1) and kbox_j5.json (lam5):
   roots equal kbox_run.roots(j, nt) (they tile B_j), every split is a midpoint
   bisection, every leaf has kcert.F_lower >= 0.
 5 ledger.py arithmetic with the certified lam1, lam5; U_new ceiling < 0.629526945573.
Run from any folder (stage5_kernelbox/run.py calls it):
  nice -n 10 python3 src/check_kernelbox.py --workers 4
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse, json, subprocess, sys, time
import multiprocessing as mp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import paths
import kcert as kc
import kformula as kf
import kbox_run as kr

AA = Q(532, 1000)


def gate1():
    a2 = AA * AA
    c4, c2, c0 = Q(1, 16), Q(3, 16) - a2 / 4, Q(1, 8) - a2 / 4
    assert c4 >= 0 and c2 >= 0 and c0 >= 0
    # polynomial identity of degree 4: equality at 5+ points suffices; use 9
    for t in [Q(k, 4) for k in range(-4, 5)]:
        assert t * t * (t * t - a2) == c4 * (16 * t ** 4 - 12 * t * t + 1) + c2 * (4 * t * t - 1) + c0
    for R in kf.BREAKS[1:6]:
        assert 1 - 1 / (2 * R * R) <= AA
    assert 2 + 2 * AA > Q(49, 16)
    return dict(c4=str(c4), c2=str(c2), c0=str(c0))


def gate3():
    pts = [(Q(1), Q(1), Q(-864, 1000)), (Q(1), kf.BREAKS[1], Q(-7, 10)), (kf.BREAKS[1], kf.BREAKS[1], Q(-99, 100)),
           (Q(1), kf.BREAKS[5], Q(-6, 10)), ((kf.BREAKS[2] + kf.BREAKS[3]) / 2, kf.BREAKS[4], Q(-8, 10)),
           (kf.BREAKS[5], kf.BREAKS[5], Q(-55, 100))]
    out = []
    for r, s, t in pts:
        lb = kc.F_lower(((r, r), (s, s), (t, t)), Q(0), AA)
        ref = kf.minusK_point(r, s, t)
        # a point box: the lower bound must not exceed -K and must agree to 1e-30
        assert lb <= ref.upper() and ref.lower() - lb < kf.A(Q(1, 10 ** 30)), (r, s, t)
        out.append([str(r), str(s), str(t), float(ref.mid())])
    return out


def replay(args):
    enc, lam = args
    br = tuple(map(Q, enc["r"])); bs = tuple(map(Q, enc["s"])); bt = tuple(map(Q, enc["t"]))
    open_r, open_s = enc["open_r"], enc["open_s"]
    stats = dict(leaves=0, min_lower=None, bad=None)

    def walk(box, node):
        if stats["bad"]:
            return
        if node == "L":
            lo = kc.F_lower(box, lam, AA, open_r, open_s)
            if not (lo >= 0):
                stats["bad"] = ["negative leaf", [[str(x) for x in iv] for iv in box]]
                return
            stats["leaves"] += 1
            v = float(lo.mid())
            stats["min_lower"] = v if stats["min_lower"] is None else min(stats["min_lower"], v)
            return
        if not (isinstance(node, list) and len(node) == 3 and node[0] in (0, 1, 2)):
            stats["bad"] = ["bad node", str(node)[:80]]
            return
        d = node[0]; a, b = box[d]; m = (a + b) / 2
        left = list(box); right = list(box); left[d] = (a, m); right[d] = (m, b)
        walk(tuple(left), node[1]); walk(tuple(right), node[2])

    walk((br, bs, bt), enc["tree"])
    return stats


def gate4(path, j, lam, workers):
    data = json.loads(Path(path).read_text())
    summ = data["summary"]
    assert summ["j"] == j and Q(summ["lam"]) == lam and summ["status"] == "ALL_LEAVES_NONNEGATIVE"
    expect = kr.roots(j, summ["nt"])
    got = [((Q(e["r"][0]), Q(e["r"][1])), (Q(e["s"][0]), Q(e["s"][1])), (Q(e["t"][0]), Q(e["t"][1])),
            e["open_r"], e["open_s"]) for e in data["roots"]]
    assert got == expect, "roots do not match the tiling of B_j"
    # tiling facts: t chunks cover [-1,-a]; regions cover [1, R_j]
    assert expect[0][2][0] == -1 and max(x[2][1] for x in expect) == -AA
    assert {x[0] for x in expect} == {(kf.BREAKS[k], kf.BREAKS[k + 1]) for k in range(j)}
    with mp.Pool(workers) as pool:
        res = pool.map(replay, [(e, lam) for e in data["roots"]], chunksize=4)
    bad = [r["bad"] for r in res if r["bad"]]
    assert not bad, bad[:3]
    return dict(file=str(path), j=j, lam=str(lam), roots=len(res), leaves=sum(r["leaves"] for r in res),
                min_leaf_lower=min(r["min_lower"] for r in res))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--lam1", default="3/1000000")
    ap.add_argument("--lam5", default="1/1000000")
    ap.add_argument("--j1-file", default=str(paths.KBOX / "kbox_j1.json"))
    ap.add_argument("--j5-file", default=str(paths.KBOX / "kbox_j5.json"))
    ap.add_argument("--out", default=str(paths.OUT / "check_kernelbox.json"), help="controls write elsewhere")
    args = ap.parse_args()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time(); out = dict(status="RUNNING")
    try:
        out["gate1"] = gate1()
        r = subprocess.run([sys.executable, str(HERE / "validate_kformula.py")], capture_output=True, text=True)
        assert r.returncode == 0 and "KFORMULA_MATCHES_RAW_V7 72" in r.stdout, r.stdout + r.stderr
        out["gate2"] = r.stdout.strip()
        out["gate3"] = gate3()
        out["gate4_j1"] = gate4(args.j1_file, 1, Q(args.lam1), args.workers)
        out["gate4_j5"] = gate4(args.j5_file, 5, Q(args.lam5), args.workers)
        r = subprocess.run([sys.executable, str(HERE / "ledger.py"), "--lam1", args.lam1, "--lam5", args.lam5],
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        led = json.loads(r.stdout)
        assert Q(led["U_new_strict_ceiling"]) < Q(629526945573, 10 ** 12)
        out["gate5"] = {k: led[k] for k in ("a_new_decimal", "a_gain", "U_new_exact_decimal",
                                            "U_new_strict_ceiling", "U_new_ceiling_decimal", "improvement")}
        out["gate5"]["binding"] = led["binding"][:2]
        out["status"] = "KERNELBOX_CERTIFIED"
    except Exception as e:  # fail closed
        out["status"] = "FAILED"; out["error"] = repr(e)[:2000]
    out["seconds"] = time.time() - t0
    Path(args.out).write_text(json.dumps(out, indent=1) + "\n")
    print(out["status"], f"{out['seconds']:.0f}s")
    sys.exit(0 if out["status"] == "KERNELBOX_CERTIFIED" else 1)
