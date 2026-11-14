"""Run complete registered wrappers; preserve a cumulative campaign ledger."""

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse, json, subprocess, sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("dimensions", type=int, nargs="+")
ap.add_argument("--jobs", type=int, default=2)
args = ap.parse_args()
report = HERE / "REGISTERED_REPLAY.json"
checks = json.loads(report.read_text())["checks"] if report.exists() else []
checks = [r for r in checks if r["dimension"] not in args.dimensions]


def run(n):
    D = ROOT / f"verifications/dim{n}/refinement_combined20261115"
    result = subprocess.run(
        [sys.executable, str(D / "verify.py")], capture_output=True, text=True
    )
    (D / "registered_replay.log").write_text(result.stdout + result.stderr)
    assert result.returncode == 0, (n, result.stderr)
    print("REGISTERED PASS", n, flush=True)
    return {"dimension": n, "registered_wrapper_pass": True}


with ThreadPoolExecutor(max_workers=args.jobs) as pool:
    for future in as_completed([pool.submit(run, n) for n in args.dimensions]):
        checks.append(future.result())
report.write_text(
    json.dumps(
        {"all_passed": True, "checks": sorted(checks, key=lambda r: r["dimension"])},
        indent=2,
    )
    + "\n"
)
