#!/usr/bin/env python3
"""Run every check behind the D4 bound and compose the final inequality.

    python3 verify/verify_all.py --out DIR [--mode quick|full] [--workers N]

Run from the D4 folder (the folder that holds paper.tex, data/ and verify/).
DIR must not exist. Each stage writes DIR/<stage>/RESULT.json; this driver passes each
stage's fresh output to the next stage, then recomputes the final ledger in exact
arithmetic from the fresh stage outputs and checks

    U_raw - a_new/8 - 3/8000000 < 629521681653/10^12.

Modes:
  full   reruns everything, including all 112 cyclic triangle branches of stage 1.
  quick  skips the 112 branch subdivisions of stage 1 (it still reruns every exact check
         and the pair pieces); stage 2 then uses the committed stage 1 values.
Two exact kernel identity checks also run in both modes.  The Lean formalization is
unreleased and is deliberately not imported by this driver.

Exit 0 means every stage passed and the final inequality holds. Anything else is failure.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import json
import os
import subprocess
import sys
import time
import hashlib

D4 = Path(__file__).resolve().parents[1]
V = D4 / "verify"
THRESHOLD = Q(629521681653, 10**12)
GRAPH_SAVING = Q(3, 8000000)


def exact(value):
    """Reject JSON floats and booleans before any rational conversion."""
    if type(value) not in (str, int):
        raise ValueError(f"expected an exact rational string or integer, got {value!r}")
    return Q(value)


def integer(value, lo, hi):
    q = exact(value)
    if q.denominator != 1 or not lo <= q <= hi:
        raise ValueError(f"expected integer in [{lo}, {hi}], got {value!r}")
    return q.numerator


def census_groups(profile):
    expected = {(0, 1)} | {(n, j) for n in range(1, 25) for j in range(1, 6)} | {(25, 4), (25, 5), (26, 5)}
    groups = {}
    for g in profile:
        key = (integer(g["N"], 0, 26), integer(g["first_radius"], 1, 5))
        if key in groups:
            raise ValueError(f"duplicate census group {key}")
        groups[key] = exact(g["new_floor"])
        if groups[key] < 0:
            raise ValueError(f"negative census floor at {key}")
    if set(groups) != expected:
        raise ValueError(f"census group domain differs: missing={expected-set(groups)}, extra={set(groups)-expected}")
    return groups


def check_export_manifests():
    """Fail closed if any shipped manifest entry is missing or stale."""
    checked = 0
    for manifest_path in sorted(V.rglob("*MANIFEST.json")):
        manifest = json.loads(manifest_path.read_text())
        rows = manifest.get("files", [])
        if isinstance(rows, dict):
            rows = [dict(value, path=path) for path, value in rows.items()]
        for row in rows:
            relative = Path(row["path"])
            path = manifest_path.parent / relative
            if relative.is_absolute() or not path.is_file() or not path.resolve().is_relative_to(V.resolve()):
                raise ValueError(f"missing or escaping manifest entry: {manifest_path}:{relative}")
            data = path.read_bytes()
            if len(data) != row["bytes"] or hashlib.sha256(data).hexdigest() != row["sha256"]:
                raise ValueError(f"manifest mismatch: {manifest_path}:{relative}")
            checked += 1
    return checked


def fail(msg):
    print("FAIL:", msg, flush=True)
    sys.exit(1)


def run_stage(name, argv, out, log_dir):
    """Run one stage; require exit 0 and RESULT.json with the right stage name and PASS."""
    log = log_dir / f"{name}.log"
    print(f"[{time.strftime('%H:%M:%S')}] {name}: {' '.join(argv[1:])}", flush=True)
    t0 = time.time()
    with log.open("w") as fh:
        rc = subprocess.run(argv, cwd=D4, stdout=fh, stderr=subprocess.STDOUT).returncode
    secs = time.time() - t0
    res_file = out / "RESULT.json"
    if not res_file.exists():
        fail(f"{name}: no RESULT.json (exit {rc}); see {log}")
    res = json.loads(res_file.read_text())
    if rc != 0 or res.get("status") != "PASS" or res.get("stage") != name:
        fail(f"{name}: exit {rc}, status {res.get('status')}, stage {res.get('stage')}; "
             f"error: {res.get('error')}; see {log}")
    print(f"    PASS in {secs:.0f} s", flush=True)
    return res


def lam(j):
    return Q(3, 10**6) if j == 1 else Q(1, 10**6)


def QN(n):
    return Q(n * (13561 * n - 179244), 250000)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, required=True, help="new output folder")
    ap.add_argument("--mode", choices=["quick", "full"], default="full")
    ap.add_argument("--workers", type=int, choices=range(1, 5), default=2, help="workers for stages 1 and 5 (2 needs about 1 GB)")
    args = ap.parse_args()
    if sys.version_info < (3, 12):
        ap.error("Python 3.12 or later is required because NumPy 2.5.2 is pinned")
    if not __debug__ or os.environ.get("PYTHONOPTIMIZE"):
        ap.error("optimized Python disables assertions in inherited checkers; unset PYTHONOPTIMIZE and omit -O")
    out = args.out.expanduser().resolve()
    if out.exists():
        ap.error("--out must not exist")
    out.mkdir(parents=True)
    logs = out / "logs"
    logs.mkdir()
    try:
        manifest_files = check_export_manifests()
    except (KeyError, OSError, TypeError, ValueError) as exc:
        fail(f"export manifest audit failed: {exc}")
    py = sys.executable
    w = str(args.workers)
    d = {s: out / s for s in ["stage1_v7", "stage2_origin", "stage3_census",
                              "stage4_graph_slack", "stage5_kernelbox",
                              "kernel_independent", "kernel_link"]}
    t0 = time.time()

    s1_argv = [py, str(V / "stage1_v7/run.py"), "--out", str(d["stage1_v7"]), "--workers", w]
    if args.mode == "quick":
        s1_argv.append("--quick")
    s1 = run_stage("stage1_v7", s1_argv, d["stage1_v7"], logs)
    if args.mode == "full" and (s1.get("mode") != "full" or s1.get("branches_rerun") != 112):
        fail("stage1_v7 did not rerun all 112 branches in full mode")

    s2_argv = [py, str(V / "stage2_origin/run.py"), "--out", str(d["stage2_origin"])]
    if args.mode == "full":
        s2_argv += ["--stage1", str(d["stage1_v7"])]
    s2 = run_stage("stage2_origin", s2_argv, d["stage2_origin"], logs)
    s3 = run_stage("stage3_census", [py, str(V / "stage3_census/run.py"), "--out", str(d["stage3_census"]),
                                     "--stage2", str(d["stage2_origin"])], d["stage3_census"], logs)
    up = ["--census", str(d["stage3_census"]), "--stage2", str(d["stage2_origin"])]
    s4 = run_stage("stage4_graph_slack", [py, str(V / "stage4_graph_slack/run.py"),
                                          "--out", str(d["stage4_graph_slack"]), *up], d["stage4_graph_slack"], logs)
    s5 = run_stage("stage5_kernelbox", [py, str(V / "stage5_kernelbox/run.py"),
                                        "--out", str(d["stage5_kernelbox"]), *up, "--workers", w],
                   d["stage5_kernelbox"], logs)

    ki = run_stage("kernel_independent", [py, str(V / "kernel_independent/run.py"),
                   "--out", str(d["kernel_independent"]), "--workers", "1"],
                   d["kernel_independent"], logs)
    counts = ki["counts"]
    if (counts.get("full_equal"), counts.get("outside_equal"), counts.get("pair_pieces_equal")) != (84, 28, 7):
        fail("independent kernel comparison did not cover all 112 branches and 7 pair pieces")
    if ki.get("component_contracts", {}).get("status") != "PASS":
        fail("independent kernel component contracts failed")
    kl = run_stage("kernel_link", [py, str(V / "kernel_link/run.py"), "--out", str(d["kernel_link"])],
                   d["kernel_link"], logs)
    if (kl.get("ordered_radial_cells"), kl.get("endpoint_cells"), kl.get("stage1_pair_pieces")) != (25, 25, 7):
        fail("stage 1/stage 5 kernel link lacks complete radial/endpoint coverage")
    if kl.get("includes_B1") is not True or kl.get("includes_B5") is not True:
        fail("kernel link does not cover both certificate domains")

    # Independent composition from the fresh outputs.
    u_raw = {"stage1 raw_objective": exact(s1["raw_objective"]),
             "stage2 exact_source_objective": exact(s2["exact_source_objective"]),
             "stage3 exact_raw_objective": exact(s3["exact_raw_objective"]),
             "stage5 U_raw": exact(s5["U_raw"]),
             "independent U_raw": exact(ki["U_raw"])}
    if len(set(u_raw.values())) != 1:
        fail(f"U_raw disagrees between stages: {u_raw}")
    U_raw = u_raw["stage1 raw_objective"]
    if args.mode == "full" and s2.get("stage1_checked") is not True:
        fail("stage2_origin did not check the fresh stage 1 output")
    if s3.get("stage2_fresh") is not True:
        fail("stage3_census did not use the fresh stage 2 output")
    if Q(s3["graph_saving"]) != GRAPH_SAVING or Q(s4["density_saving"]) != GRAPH_SAVING:
        fail("graph saving is not 3/8000000 in stages 3 and 4")

    census = json.loads((d["stage3_census"] / "FINITE_COUNT_ACCEPTANCE.json").read_text())
    groups = census_groups(census["profile"])
    if min(groups.values()) != Q(s3["new_beta"]):
        fail("census minimum differs from stage 3 new_beta")
    raised = {k: f + lam(k[1]) * max(QN(k[0]), Q(0)) for k, f in groups.items()}
    a_new = min(raised.values())
    binding = min(raised, key=raised.get)
    if a_new != Q(s5["a_new"]) or list(binding) != list(s5["binding_group"]):
        fail("recomputed a_new or binding group differs from stage 5")
    U_new = U_raw - a_new / 8 - GRAPH_SAVING
    if U_new != Q(s5["U_new"]):
        fail("recomputed U_new differs from stage 5")
    if not U_new < THRESHOLD:
        fail(f"U_new = {float(U_new)} is not below {THRESHOLD}")

    summary = {
        "status": "PASS", "mode": args.mode, "full_certificate_replay": args.mode == "full",
        "stages": {k: "PASS" for k in d},
        "branches_rerun": s1.get("branches_rerun"),
        "count_vectors_checked": s2.get("count_vectors_checked"),
        "group_floors": len(groups),
        "independent_kernel": {"full_branches": counts["full_equal"], "exterior_branches": counts["outside_equal"],
                               "pair_pieces": counts["pair_pieces_equal"],
                               "polynomial_coefficients": counts["polynomial_terms_compared"],
                               "nonpolynomial_identity": "written component proof with exact parameter contracts"},
        "kernel_link": {"radial_cells": kl["ordered_radial_cells"], "angular_endpoint_cells": kl["endpoint_cells"],
                        "method": kl["identity_method"], "negative_controls_rejected": kl["negative_controls_rejected"]},
        "kernel_box_leaves_replayed": [s5["kernel_box"].get("j1_leaves_replayed"),
                                       s5["kernel_box"].get("j5_leaves_replayed")],
        "U_raw": str(U_raw), "a_new": str(a_new), "binding_group": list(binding),
        "graph_saving": str(GRAPH_SAVING), "U_new": str(U_new),
        "U_new_decimal": f"{float(U_new):.15f}",
        "bound": f"density < {THRESHOLD} = 0.629521681653",
        "seconds": round(time.time() - t0, 1),
        "manifest_files_checked": manifest_files,
    }

    (out / "SUMMARY.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    if args.mode == "full":
        print("PASS: every sphere packing in R^4 has density < 0.629521681653 (full mode)")
    else:
        print("PASS: quick checks only; the 112 cyclic sign subdivisions were NOT rerun")


if __name__ == "__main__":
    main()
