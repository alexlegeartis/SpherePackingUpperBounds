"""Negative controls for the frozen wrappers' actual registry/hash prefix.

These controls do not run or replace mathematical proof programs. The complete
registered replays are performed separately.
"""

from pathlib import Path
from fractions import Fraction
import argparse, ast, json, shutil, subprocess, sys, tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("dimensions", nargs="+", type=int)
args = ap.parse_args()
runner = """from pathlib import Path
import ast,sys
file=Path(sys.argv[1]);text=file.read_text();tree=ast.parse(text)
cut=next(node.lineno for node in tree.body if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=="executed" for t in node.targets))
prefix="\\n".join(text.splitlines()[:cut-1])
exec(compile(prefix,str(file),"exec"),{"__file__":str(file),"__name__":"__main__"})
"""
report_file = HERE / "REGISTRATION_CONTROLS.json"
checks = json.loads(report_file.read_text())["checks"] if report_file.exists() else []
checks = [row for row in checks if row["dimension"] not in args.dimensions]
for n in args.dimensions:
    source = ROOT / f"verifications/dim{n}/refinement_combined20261115"
    cfg = json.loads((source / "config.json").read_text())
    for name in [
        "valid_current",
        "tighter_current",
        "weaker_current",
        "wrong_package",
        "tampered_witness",
        "tampered_core",
        "valid_historical",
    ]:
        with tempfile.TemporaryDirectory(
            prefix="packing-registration-control-"
        ) as directory:
            root = Path(directory)
            D = root / f"verifications/dim{n}/refinement_combined20261115"
            shutil.copytree(source, D)
            shared = root / "verifications/profile_enrichment_20261114"
            shutil.copytree(ROOT / "verifications/profile_enrichment_20261114", shared)
            row = {
                "dimension": n,
                "package": cfg["package"],
                "strict_density_upper_bound": cfg["strict_density_upper_bound"],
            }
            expected = name in ["valid_current", "valid_historical"]
            if name == "tighter_current":
                row["strict_density_upper_bound"] = str(
                    Fraction(cfg["strict_density_upper_bound"]) / 2
                )
            elif name == "weaker_current":
                row["strict_density_upper_bound"] = str(
                    Fraction(cfg["strict_density_upper_bound"]) * 2
                )
            elif name == "wrong_package":
                row["package"] = "verifications/unverified"
            elif name == "tampered_witness":
                (D / "certificate.json").write_text(
                    (D / "certificate.json").read_text() + " "
                )
            elif name == "tampered_core":
                (shared / "core.py").write_text(
                    (shared / "core.py").read_text() + "\n# altered\n"
                )
            elif name == "valid_historical":
                row.update(
                    package="verifications/a_separately_registered_successor",
                    earlier_repository_packages=[cfg["package"]],
                    strict_density_upper_bound=str(
                        Fraction(cfg["strict_density_upper_bound"]) / 2
                    ),
                )
            (root / "RESULTS.json").write_text(json.dumps([row]))
            result = subprocess.run(
                [sys.executable, "-c", runner, str(D / "verify.py")],
                capture_output=True,
                text=True,
            )
            assert (result.returncode == 0) == expected, (
                n,
                name,
                result.stdout,
                result.stderr,
            )
            checks.append(
                {
                    "dimension": n,
                    "control": name,
                    "expected_prefix_pass": expected,
                    "observed_prefix_pass": result.returncode == 0,
                }
            )
(HERE / "REGISTRATION_CONTROLS.json").write_text(
    json.dumps(
        {
            "all_passed": True,
            "scope": "Actual source-hash/registry prefix only; mathematical proofs are replayed separately.",
            "checks": checks,
        },
        indent=2,
    )
    + "\n"
)
print("PASS:", len(checks), "registry/source-hash controls")
