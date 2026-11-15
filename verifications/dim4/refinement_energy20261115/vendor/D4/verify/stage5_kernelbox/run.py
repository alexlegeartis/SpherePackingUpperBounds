"""Stage 5: far-pair kernel box and final ledger. Fail closed.

Runs, in order, each as a child process (any nonzero exit, timeout or missing
output is FAIL):
  1 src/check_kernelbox.py  gates 1-5: exact identities, kformula against the
    independent raw_v7 evaluator at 72 points, point consistency, replay of
    every leaf of data/kernel-box/kbox_j1.json and kbox_j5.json, ledger.
  2 src/hrlemma_check.py    exact facts behind H_R >= lambda_j max(Q_N, 0).
  3 src/coord_check_kernelbox.py  raw_v7 additivity on 3- and 4-point stars and
    -K - lambda_1 q >= 0 at 88 rational points of B_1 (independent of kformula).
  4 src/ledger.py           a_new, U_new and the strict ceiling.
  5 src/exact_values.py     exact a_new, U_new, next floor (26,5).
Then this script recomputes a_new and U_new itself from the census floors and
compares everything with data/kernel-box/ledger.json.

Usage: python3 run.py --out DIR [--workers 4] [--census DIR] [--stage2 DIR]
--census / --stage2 take fresh upstream output folders; their RESULT.json must say PASS.
Writes DIR/RESULT.json; exit 0 iff status PASS.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse, hashlib, json, os, subprocess, sys, threading, time

sys.set_int_max_str_digits(0)
STAGE = Path(__file__).resolve().parent
SRC = STAGE / "src"
PKG = STAGE.parents[1]
DATA = PKG / "data"
LAM1, LAM5 = Q(3, 10 ** 6), Q(1, 10 ** 6)
A = Q(133, 250)
GRAPH_SAVING = Q(3, 8000000)
CEILING = Q(629521681653, 10 ** 12)
TIMEOUT = 3600


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


class RSS:
    """Samples the summed RSS of this process and all its descendants every 2 s."""

    def __init__(self):
        self.peak = 0; self.stop = False
        self.t = threading.Thread(target=self.loop, daemon=True); self.t.start()

    def loop(self):
        me = os.getpid()
        while not self.stop:
            try:
                out = subprocess.run(["ps", "-A", "-o", "pid=,ppid=,rss="], capture_output=True, text=True).stdout
                rows = [tuple(map(int, l.split())) for l in out.splitlines() if l.strip()]
                kids = {}
                for pid, ppid, rss in rows:
                    kids.setdefault(ppid, []).append((pid, rss))
                own = {pid: rss for pid, _, rss in rows}
                total, stack = own.get(me, 0), [me]
                while stack:
                    for pid, rss in kids.get(stack.pop(), []):
                        total += rss; stack.append(pid)
                self.peak = max(self.peak, total)
            except Exception:
                pass
            time.sleep(2)


def child(name, args, env, log, detail=None):
    """Run one script. detail: its JSON output, whose "error" is copied into a failure message."""
    t0 = time.time()
    try:
        r = subprocess.run([sys.executable] + args, capture_output=True, text=True, env=env, timeout=TIMEOUT, cwd=STAGE)
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"{name}: timeout after {TIMEOUT} s")
    log[name] = dict(exit=r.returncode, seconds=round(time.time() - t0, 1), stdout_tail=r.stdout[-600:], stderr_tail=r.stderr[-600:])
    if r.returncode != 0:
        why = ""
        if detail is not None and Path(detail).exists():
            try:
                why = json.loads(Path(detail).read_text()).get("error", "")
            except Exception:
                pass
        raise RuntimeError(f"{name}: exit {r.returncode}: {(r.stdout + r.stderr)[-800:]} {why}")
    return r.stdout


def upstream(census, stage2):
    """Resolve the census and stage-2 files. A folder given on the command line is a fresh
    upstream output: its RESULT.json must say PASS for the right stage, or this stage fails."""
    def passed(d, stage):
        r = d / "RESULT.json"
        if not r.exists():
            raise RuntimeError(f"missing {r}")
        res = json.loads(r.read_text())
        if res.get("stage") != stage or res.get("status") != "PASS":
            raise RuntimeError(f"{r}: stage {res.get('stage')!r} status {res.get('status')!r}, need {stage!r} PASS")
    if census:
        d = Path(census).resolve(); passed(d, "stage3_census")
        names = [n for n in ("census_floors.json", "FINITE_COUNT_ACCEPTANCE.json") if (d / n).exists()]
        if not names:
            raise RuntimeError(f"no census_floors.json or FINITE_COUNT_ACCEPTANCE.json in {d}")
        census_file = d / names[0]
    else:
        census_file = DATA / "census_floors.json"
    if stage2:
        d = Path(stage2).resolve(); passed(d, "stage2_origin")
    else:
        d = PKG / "verify/stage2_origin/reference"
    return census_file, d / "ACCEPTANCE.json"


def load(p):
    p = Path(p)
    if not p.exists():
        raise RuntimeError(f"missing output {p}")
    return json.loads(p.read_text())


def QN(N):
    return Q(N * (13561 * N - 179244), 250000)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--census", help="fresh stage 3 output folder: RESULT.json PASS plus census_floors.json or FINITE_COUNT_ACCEPTANCE.json (default: data/census_floors.json)")
    ap.add_argument("--stage2", help="fresh stage 2 output folder: RESULT.json PASS plus ACCEPTANCE.json (default: verify/stage2_origin/reference)")
    ap.add_argument("--lam1", default=str(LAM1), help="negative controls only")
    ap.add_argument("--lam5", default=str(LAM5), help="negative controls only")
    args = ap.parse_args()
    out = Path(args.out).resolve(); out.mkdir(parents=True, exist_ok=True)
    lam1, lam5 = Q(args.lam1), Q(args.lam5)
    res = dict(stage="stage5_kernelbox", status="FAIL", lam1=str(lam1), lam5=str(lam5), a=str(A))
    log = {}; rss = RSS(); t0 = time.time()
    try:
        census_file, acc_file = upstream(args.census, args.stage2)
        stage2_dir = acc_file.parent
        res["census_file"] = str(census_file); res["stage2_file"] = str(acc_file)
        env = dict(os.environ, D4_STAGE5_OUT=str(out), D4_STAGE2_DIR=str(stage2_dir),
                   D4_CENSUS_FILE=str(census_file), PYTHONDONTWRITEBYTECODE="1")
        import flint
        res["python"] = sys.version.split()[0]; res["python_flint"] = flint.__version__
        inputs = [census_file, stage2_dir / "ACCEPTANCE.json", DATA / "kernel-box/kbox_j1.json",
                  DATA / "kernel-box/kbox_j5.json", DATA / "kernel-box/ledger.json"] + \
                 [DATA / "v7" / f"{n}.json" for n in ("kernel", "origin", "scalar", "graph_energies")]
        res["inputs"] = {str(p.relative_to(PKG)) if p.is_relative_to(PKG) else str(p): sha(p) for p in inputs}

        # 1 kernel box: identities, kformula vs raw_v7, replay of every leaf, ledger
        child("check_kernelbox", [str(SRC / "check_kernelbox.py"), "--workers", str(args.workers), "--lam1", args.lam1,
                                  "--lam5", args.lam5, "--out", str(out / "check_kernelbox.json")], env, log,
              detail=out / "check_kernelbox.json")
        ck = load(out / "check_kernelbox.json")
        assert ck["status"] == "KERNELBOX_CERTIFIED", ck.get("error")
        assert ck["gate2"].startswith("KFORMULA_MATCHES_RAW_V7 72"), ck["gate2"]
        vk = load(out / "validate_kformula.json")
        assert vk["status"] == "KFORMULA_MATCHES_RAW_V7" and vk["points"] == 72 and all(r["ok"] for r in vk["rows"])
        assert len(vk["q0"]) == 9
        g1, g5 = ck["gate4_j1"], ck["gate4_j5"]
        assert (g1["j"], Q(g1["lam"]), g5["j"], Q(g5["lam"])) == (1, lam1, 5, lam5)
        assert g1["roots"] == 94 and g5["roots"] == 2350 and g1["leaves"] > 0 and g5["leaves"] > 0
        res["kernel_box"] = dict(j1_roots=g1["roots"], j1_leaves_replayed=g1["leaves"], j1_min_leaf_lower=g1["min_leaf_lower"],
                                 j5_roots=g5["roots"], j5_leaves_replayed=g5["leaves"], j5_min_leaf_lower=g5["min_leaf_lower"],
                                 kformula_vs_raw_v7_points=vk["points"], q0_radii=len(vk["q0"]),
                                 q_gegenbauer=[ck["gate1"]["c4"], ck["gate1"]["c2"], ck["gate1"]["c0"]])

        # 2 far-pair lemma facts
        child("hrlemma_check", [str(SRC / "hrlemma_check.py")], env, log)
        hr = load(out / "hrlemma_check.json")
        assert hr["status"] == "PASS" and all(v["ok"] for k, v in hr["checks"].items() if isinstance(v, dict) and "ok" in v)

        # 3 independent spot checks through raw_v7 only
        child("coord_check_kernelbox", [str(SRC / "coord_check_kernelbox.py")], env, log)
        cc = load(out / "coord_check_kernelbox.json")["summary"]
        assert cc["additivity_all_overlap"] is True and cc["band_points"] >= 80 and cc["band_failures"] == 0, cc
        res["raw_v7_spot_checks"] = cc

        # 4 ledger, compared with the committed ledger
        child("ledger", [str(SRC / "ledger.py"), "--lam1", args.lam1, "--lam5", args.lam5, "--out", str(out / "ledger.json")], env, log)
        led = load(out / "ledger.json")
        committed = load(DATA / "kernel-box/ledger.json")
        res["ledger_equals_committed"] = led == committed
        assert led == committed, "recomputed ledger differs from data/kernel-box/ledger.json"

        # 5 exact values
        child("exact_values", [str(SRC / "exact_values.py")], env, log)
        ev = load(out / "exact_values.json")
        assert ev["status"] == "COMPUTED_EXACT"

        # 6 independent recomputation here, from the census floors only
        census = load(census_file)
        acc = load(stage2_dir / "ACCEPTANCE.json")
        U_raw = Q(census["exact_raw_objective"])
        assert U_raw == Q(acc["exact_source_objective"]), "U_raw differs between census and stage 2"
        assert Q(census["graph_saving"]) == GRAPH_SAVING
        prof = census["profile"]
        keys = {(g["N"], g["first_radius"]) for g in prof}
        assert len(prof) == 124 and len(keys) == 124 and all(g["first_radius"] in (1, 2, 3, 4, 5) for g in prof)
        raised = sorted((Q(g["new_floor"]) + (lam1 if g["first_radius"] == 1 else lam5) * max(QN(g["N"]), Q(0)), g["N"], g["first_radius"])
                        for g in prof)
        a_new, N0, j0 = raised[0]
        U_new = U_raw - a_new / 8 - GRAPH_SAVING
        assert a_new == Q(led["a_new"]) == Q(ev["a_new"])
        assert U_new == Q(ev["U_new"])
        assert U_new < CEILING, "U_new is not below 0.629521681653"
        grid = Q(-((-U_new.numerator * 10 ** 12) // U_new.denominator), 10 ** 12)
        assert grid == CEILING == Q(led["U_new_strict_ceiling"])
        res.update(binding_group=[N0, j0], Q_binding=str(QN(N0)), a_plus=census["new_beta"], a_new=str(a_new),
                   U_raw=str(U_raw), graph_saving=str(GRAPH_SAVING), U_new=str(U_new),
                   U_new_decimal=f"{float(U_new):.15f}", ceiling=str(CEILING), U_new_below_ceiling=True,
                   next_group=[raised[1][1], raised[1][2]], next_raised_floor=str(raised[1][0]))
        res["status"] = "PASS"
    except Exception as e:  # fail closed
        res["status"] = "FAIL"; res["error"] = repr(e)[:3000]
    rss.stop = True
    res["children"] = log
    res["seconds"] = round(time.time() - t0, 1)
    res["peak_total_rss_mb"] = round(rss.peak / 1024, 1)
    (out / "RESULT.json").write_text(json.dumps(res, indent=1) + "\n")
    print(res["status"], res.get("error", ""), f"{res['seconds']} s")
    sys.exit(0 if res["status"] == "PASS" else 1)


if __name__ == "__main__":
    main()
