"""Run all four frozen proof paths before installing a registry entry."""

import argparse, concurrent.futures, json, subprocess, sys, hashlib, time
from pathlib import Path

if not __debug__:
    raise RuntimeError("Assertions must remain enabled")
ROOT = Path(__file__).resolve().parents[2]
ap = argparse.ArgumentParser()
ap.add_argument("dimensions", type=int, nargs="+")
ap.add_argument("--jobs", type=int, default=3)
args = ap.parse_args()
assert 1 <= args.jobs <= 4


def replay(n):
    D = ROOT / f"verifications/dim{n}/refinement_combined20261115"
    manifest = json.loads((D / "source_hashes.json").read_text())
    for name, digest in manifest.items():
        assert hashlib.sha256((D / name).read_bytes()).hexdigest() == digest
    start = time.time()
    executed = []
    for script in [
        "primary_angular.py",
        "primary_radial.py",
        "independent_angular.py",
        "independent_radial.py",
    ]:
        result = subprocess.run(
            [
                sys.executable,
                str(ROOT / "verifications/profile_enrichment_20261114" / script),
                str(D),
            ],
            capture_output=True,
            text=True,
        )
        (D / (script.removesuffix(".py") + ".log")).write_text(
            result.stdout + result.stderr
        )
        print(n, script, result.returncode, result.stdout[-240:], flush=True)
        assert result.returncode == 0, (n, script, result.stderr)
        executed.append({"script": script, "returncode": 0})
    cfg = json.loads((D / "config.json").read_text())
    out = {
        "pass": True,
        "dimension": n,
        "strict_density_upper_bound": cfg["strict_density_upper_bound"],
        "complete_primary_and_independent_replay": True,
        "executed": executed,
        "source_manifest_sha256": hashlib.sha256(
            (D / "source_hashes.json").read_bytes()
        ).hexdigest(),
        "external_neighbor_counts_required": False,
        "seconds": time.time() - start,
        "arithmetic_backend": "FLINT/Arb shared; independent proof implementations",
        "execution": "All four scripts invoked directly before registry installation; wrapper registry/hash gates are checked separately",
    }
    (D / "full_replay_result.json").write_text(json.dumps(out, indent=2) + "\n")
    print("FOUR PATHS PASS", n, cfg["strict_density_upper_bound"], flush=True)


with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
    list(pool.map(replay, args.dimensions))
