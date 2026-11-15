"""Negative controls for stage 5. Each control copies the whole D4 folder,
breaks one thing in the copy, runs run.py there and passes only if run.py
fails closed (exit 1, RESULT.json status FAIL) for the expected reason.

  lam1_4e-6        kbox_j1.json replayed at lambda_1 = 4/10^6 (summary lam
                   set to match, run.py --lam1 4/1000000). False at the dip.
  merge_j5_leaves  in kbox_j5.json, two sibling leaves of the dip root are
                   merged into their parent box (one leaf replaces a split).

Usage: python3 negative_controls.py --work DIR [--workers 2]
Writes DIR/controls.json; exit 0 iff every control fails closed as expected.
"""
from pathlib import Path
import argparse, json, shutil, subprocess, sys, time

STAGE = Path(__file__).resolve().parent
PKG = STAGE.parents[1]
DIP_T = -0.864  # t of the binding dip at radii (1, 1)


def copy_pkg(dest):
    # Proof sources and certificates suffice; a Lean dependency cache can be many GB.
    shutil.copytree(PKG, dest, ignore=shutil.ignore_patterns(
        "__pycache__", ".lake", ".git", ".DS_Store", "out"))
    return dest


def dip_root(roots):
    """Index of the root with both radii in [1, R_1] whose t-chunk contains DIP_T."""
    from fractions import Fraction as Q
    for i, e in enumerate(roots):
        if Q(e["r"][0]) == 1 and Q(e["s"][0]) == 1 and Q(e["t"][0]) <= Q(DIP_T) <= Q(e["t"][1]):
            return i
    raise SystemExit("dip root not found")


def merge_deepest(node):
    """Replace one deepest [d, 'L', 'L'] by 'L'. Returns (new_node, depth) or (node, None)."""
    best = [None, -1]

    def find(n, path, depth):
        if isinstance(n, list):
            if n[1] == "L" and n[2] == "L" and depth > best[1]:
                best[0], best[1] = list(path), depth
            find(n[1], path + [1], depth + 1); find(n[2], path + [2], depth + 1)
    find(node, [], 0)
    if best[0] is None:
        return node, None
    if not best[0]:
        return "L", 0
    parent = node
    for k in best[0][:-1]:
        parent = parent[k]
    parent[best[0][-1]] = "L"
    return node, best[1]


def run(copy, out, extra, workers):
    t0 = time.time()
    r = subprocess.run([sys.executable, str(copy / "verify/stage5_kernelbox/run.py"), "--out", str(out),
                        "--workers", str(workers)] + extra, capture_output=True, text=True)
    res = json.loads((out / "RESULT.json").read_text()) if (out / "RESULT.json").exists() else {}
    return r.returncode, res, round(time.time() - t0, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    ap.add_argument("--workers", type=int, choices=range(1, 5), default=2)
    args = ap.parse_args()
    work = Path(args.work).resolve(); work.mkdir(parents=True, exist_ok=False)
    report = {}

    # control 1: lambda_1 = 4/10^6
    c = copy_pkg(work / "lam1_4e-6" / "D4")
    p = c / "data/kernel-box/kbox_j1.json"
    d = json.loads(p.read_text()); d["summary"]["lam"] = "4/1000000"; p.write_text(json.dumps(d))
    code, res, sec = run(c, work / "lam1_4e-6" / "out", ["--lam1", "4/1000000"], args.workers)
    err = res.get("error", "")
    report["lam1_4e-6"] = dict(exit=code, status=res.get("status"), seconds=sec, error=err[:600],
                               control_pass=code == 1 and res.get("status") == "FAIL" and "negative leaf" in err)

    # control 2: merge two sibling leaves of the dip root in kbox_j5.json
    c = copy_pkg(work / "merge_j5_leaves" / "D4")
    p = c / "data/kernel-box/kbox_j5.json"
    d = json.loads(p.read_text()); i = dip_root(d["roots"])
    d["roots"][i]["tree"], depth = merge_deepest(d["roots"][i]["tree"])
    assert depth is not None, "dip root of kbox_j5 has no split to merge"
    p.write_text(json.dumps(d))
    code, res, sec = run(c, work / "merge_j5_leaves" / "out", [], args.workers)
    err = res.get("error", "")
    report["merge_j5_leaves"] = dict(root_index=i, merged_at_depth=depth, exit=code, status=res.get("status"),
                                     seconds=sec, error=err[:600],
                                     control_pass=code == 1 and res.get("status") == "FAIL" and "negative leaf" in err)

    ok = all(v["control_pass"] for v in report.values())
    out = dict(status="ALL_CONTROLS_FAIL_CLOSED" if ok else "CONTROL_PROBLEM", controls=report)
    (work / "controls.json").write_text(json.dumps(out, indent=1) + "\n")
    for k, v in report.items():
        print(k, "PASS" if v["control_pass"] else "UNEXPECTED", v["exit"], v["error"][:160])
    print(out["status"])
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
