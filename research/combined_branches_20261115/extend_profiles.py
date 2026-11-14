"""Apply the merged profile tools to dimensions omitted by the earlier run.

Each stage runs separately. Failed numerical/exact stages cannot register a
bound. Successful exports are frozen and replayed before later registration.
"""

from pathlib import Path
from fractions import Fraction
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse, json, os, subprocess, sys, time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TOOLS = ROOT / "research/ambitious_profiles_20261114"
ap = argparse.ArgumentParser()
ap.add_argument("dimensions", type=int, nargs="+")
ap.add_argument("--jobs", type=int, default=2)
args = ap.parse_args()
logs = HERE / "logs"
logs.mkdir(exist_ok=True)
rows = json.loads((ROOT / "RESULTS.json").read_text())
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")


def run(n):
    started = time.time()
    angular = f"integrated_n{n}_a-96_b120_ang19"
    scale = Fraction(n, 8)
    tag = f'integrated_n{n}_d95_s{str(scale).replace("/","_")}_a-96_b120'
    folder = TOOLS / tag

    def command(stage, *argv):
        result = subprocess.run(
            [sys.executable, *map(str, argv)],
            cwd=ROOT,
            env=env,
            capture_output=True,
            text=True,
        )
        (logs / f"n{n}_{stage}.log").write_text(result.stdout + result.stderr)
        if result.returncode:
            raise RuntimeError(
                f"{stage}: exit {result.returncode}; see logs/n{n}_{stage}.log"
            )

    try:
        command(
            "angular",
            TOOLS / "search_angular.py",
            n,
            "--a=-96",
            "--b=120",
            "--rounds",
            6,
            "--tag",
            angular,
        )
        command(
            "angular_check",
            TOOLS / "check_angular.py",
            TOOLS / angular / "angular.json",
        )
        command(
            "radial",
            TOOLS / "search_radial.py",
            n,
            "--s",
            scale,
            "--a=-96",
            "--b=120",
            "--degree",
            95,
            "--cut",
            60,
            "--angular",
            TOOLS / angular / "angular.json",
            "--tag",
            tag,
        )
        candidate = json.loads((folder / "candidate.json").read_text())
        old = next(r for r in rows if r["dimension"] == n)
        if candidate["density_candidate"] >= float(old["strict_density_upper_bound"]):
            return {
                "dimension": n,
                "numerical_improvement": False,
                "candidate": candidate["density_candidate"],
            }
        command(
            "continue",
            TOOLS / "continue_radial.py",
            folder / "candidate.json",
            "--rounds",
            24,
        )
        command(
            "freeze", HERE / "build_packages.py", folder / "radial_certificate.json"
        )
        command("frozen_replay", HERE / "replay_frozen.py", n, "--jobs", 1)
        certificate = json.loads((folder / "radial_certificate.json").read_text())
        return {
            "dimension": n,
            "frozen_four_path_pass": True,
            "previous_strict_upper": old["strict_density_upper_bound"],
            "new_strict_upper": certificate["strict_density_upper_bound"],
            "seconds": time.time() - started,
            "discovery_directory": str(folder.relative_to(ROOT)),
        }
    except Exception as error:
        return {
            "dimension": n,
            "passed": False,
            "error": str(error),
            "seconds": time.time() - started,
        }


report_file = HERE / "profile_extension_results.json"
result = json.loads(report_file.read_text()) if report_file.exists() else []
with ThreadPoolExecutor(max_workers=args.jobs) as pool:
    futures = [pool.submit(run, n) for n in args.dimensions]
    for future in as_completed(futures):
        row = future.result()
        result = [old for old in result if old["dimension"] != row["dimension"]] + [row]
        print(json.dumps(row), flush=True)
        (HERE / "profile_extension_results.json").write_text(
            json.dumps(sorted(result, key=lambda r: r["dimension"]), indent=2) + "\n"
        )
