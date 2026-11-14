"""Replay the compact algebraic checks for this research checkpoint.

This checks the recorded exhaustive interval outputs. To repeat the complete
4.88-billion-node computation, separately run ne64_hole_certificate/run.py.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import json
import os
import subprocess
import sys
import time

P = Path(__file__).resolve().parent


def main():
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    logdir = P / "verification_logs"
    logdir.mkdir(exist_ok=True)
    checks = []

    def run(name):
        start = time.monotonic()
        logfile = logdir / (name.replace("/", "_") + ".log")
        with logfile.open("w") as out:
            result = subprocess.run([sys.executable, str(P / name)], stdout=out,
                                    stderr=subprocess.STDOUT, env=env, cwd=P)
        row = {"script": name, "exit_code": result.returncode,
               "elapsed_seconds": round(time.monotonic() - start, 3)}
        print(json.dumps(row), flush=True)
        if result.returncode:
            raise RuntimeError(f"Check failed: {name}; see {logfile}")
        return row

    # Inputs, then shard aggregation, then the dependent 65D Gram proof.
    for name in ["ne64_hole_certificate/prepare.py", "ne64_hole_certificate/summarize.py",
                 "higher_parent/verify_ne65.py"]:
        checks.append(run(name))
    independent = [
        "higher_parent/verify_gamma72.py",
        "higher_parent/ne64_theta.py",
        "higher_parent/audit_ne64_rounding.py",
        "audits/audit_ne64_antipodes.py",
        "audits/verify_ne64_sections.py",
        "audits/verify_gamma72_leech_sections.py",
        "code_routes/verify_modular_obstruction.py",
        "code_routes/verify_rank1_antipode.py",
        "code_routes/verify_ne64_theta_count.py",
        "code_routes/verify_ne64_discriminant_action.py",
        "code_routes/verify_torus_layer_bounds.py",
        "code_routes/verify_71_discriminants.py",
        "transfer_neighbors/verify_transfer.py",
        "transfer_neighbors/verify_dual_quotient.py",
    ]
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures = [pool.submit(run, name) for name in independent]
        for future in as_completed(futures):
            checks.append(future.result())
    checks.append(run("audits/propagated_benchmarks.py"))

    def load(name):
        return json.loads((P / name).read_text())

    assert load("higher_parent/ne65_verified.json")["status"] == "verified_minimum_12_lattice"
    assert load("ne64_hole_certificate/summary.json")["certified_empty"] is True
    gamma = load("higher_parent/gamma72_verified.json")
    assert gamma["all_checked"] is True and len(gamma["certificates"]) == 10
    assert set(range(66, 72)) <= {r["dimension"] for r in gamma["certificates"]}
    antipodes = load("audits/ne64_antipodes_audit.json")
    assert antipodes["all_passed"] is True and len(antipodes["checks"]) == 4
    quotient = load("transfer_neighbors/dual_quotient_verified.json")
    assert {r["index"] for r in quotient} == set(range(7))
    assert all(r["all_classes_covered"] for r in quotient)
    report = {
        "all_compact_checks_passed": True,
        "scripts_checked": len(checks),
        "checks": sorted(checks, key=lambda r: r["script"]),
        "interval_computation": "Recorded four-shard COMPLETE outputs checked; full traversal is reproduced by ne64_hole_certificate/run.py",
        "parent_minima": "Published theorems for Ne64, Gamma72, and the seven catalogue32D parents are inherited",
        "global_record_status": "Priority unestablished; comparison families are explicitly scoped",
    }
    (P / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "checks"}, indent=2))


if __name__ == "__main__":
    main()
