"""Register only fully replayed frozen refinements, preserving historical checks."""

from pathlib import Path
from fractions import Fraction
import argparse, hashlib, json, re, subprocess, sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("dimensions", type=int, nargs="+")
args = ap.parse_args()
rows = json.loads((ROOT / "RESULTS.json").read_text())
previous_report = HERE / "REGISTERED_BOUNDS.json"
changes = (
    json.loads(previous_report.read_text())["bounds"]
    if previous_report.exists()
    else []
)
assert not set(args.dimensions) & {
    row["dimension"] for row in changes
}, "Already registered in this campaign"
for n in args.dimensions:
    D = ROOT / f"verifications/dim{n}/refinement_combined20261115"
    cfg = json.loads((D / "config.json").read_text())
    replay = json.loads((D / "full_replay_result.json").read_text())
    manifest = json.loads((D / "source_hashes.json").read_text())
    assert (
        replay["pass"]
        and replay["complete_primary_and_independent_replay"]
        and replay["dimension"] == n
    )
    assert len(replay["executed"]) == 4 and all(
        x["returncode"] == 0 for x in replay["executed"]
    )
    assert replay["strict_density_upper_bound"] == cfg["strict_density_upper_bound"]
    assert (
        replay["source_manifest_sha256"]
        == hashlib.sha256((D / "source_hashes.json").read_bytes()).hexdigest()
    )
    for name, digest in manifest.items():
        assert hashlib.sha256((D / name).read_bytes()).hexdigest() == digest
    row = next(x for x in rows if x["dimension"] == n)
    old = row["strict_density_upper_bound"]
    previous_package = row["package"]
    assert Fraction(cfg["strict_density_upper_bound"]) < Fraction(old)
    earlier = row.setdefault("earlier_repository_upper_densities", [])
    if old not in earlier:
        earlier.append(old)
    packages = row.setdefault("earlier_repository_packages", [])
    if previous_package not in packages:
        packages.append(previous_package)
    row.update(
        strict_density_upper_bound=cfg["strict_density_upper_bound"],
        package=cfg["package"],
        source_package=cfg["package"],
        previous_repository_package=previous_package,
        previous_repository_strict_upper=old,
        method="Signed logarithmic radial profile and degree95 Laguerre enrichment; complete primary and independent continuum proofs",
        radial_profile={
            "a_exact": cfg["profile_a_exact"],
            "b_exact": cfg["profile_b_exact"],
            "q_exact": "1",
        },
        external_neighbor_counts_required=False,
    )
    # Historical witnesses remain independently replayable under a stronger
    # current registry bound. Their mathematical inputs and hashes are untouched.
    legacy = ROOT / f"verifications/dim{n}/verify.py"
    lines = legacy.read_text().splitlines(keepends=True)
    changed = 0
    for i, line in enumerate(lines):
        if not line.startswith("assert ") or "row['package']" not in line:
            continue
        needle = f"row['package']=='verifications/dim{n}'"
        assert needle in line, (legacy, line)
        line = line.replace(
            needle, f"row['package'] in ['verifications/dim{n}', '{cfg['package']}']"
        ).replace(
            "Fraction(row['strict_density_upper_bound'])==",
            "Fraction(row['strict_density_upper_bound'])<=",
        )
        lines[i] = line
        changed += 1
    assert changed == 1
    legacy.write_text("".join(lines))
    change = {
        "dimension": n,
        "previous_strict_upper": old,
        "strict_upper": cfg["strict_density_upper_bound"],
        "decrease_percent": float(
            100 * (1 - Fraction(cfg["strict_density_upper_bound"]) / Fraction(old))
        ),
        "package": cfg["package"],
        "four_path_pre_registration_pass": True,
    }
    changes.append(change)
    print(json.dumps(change), flush=True)
(ROOT / "RESULTS.json").write_text(json.dumps(rows, indent=2) + "\n")
runall = ROOT / "run_all.py"
text = runall.read_text()
needle = "jobs.append(('dim7','verify.py'))"
assert needle in text
text = text.replace(
    needle,
    "jobs += [(f'dim{n}', 'refinement_combined20261115/verify.py') for n in "
    + repr(sorted(args.dimensions))
    + "]\n"
    + needle,
)
runall.write_text(text)
comparison = ROOT / "comparison.py"
text = comparison.read_text()
needle = "assert set(proof_ideas)"
assert needle in text
text = text.replace(
    needle,
    "for n in "
    + repr(sorted(args.dimensions))
    + ":proof_ideas[n]='Signed logarithmic profile + degree95 radial enrichment'\n"
    + needle,
)
comparison.write_text(text)
(HERE / "REGISTERED_BOUNDS.json").write_text(
    json.dumps(
        {
            "date": "2026-11-15",
            "bounds": changes,
            "claim": "Certified upper-bound improvements over the merged main baseline; no new lower-bound record is asserted.",
        },
        indent=2,
    )
    + "\n"
)
subprocess.run(
    [sys.executable, str(ROOT / "comparison.py"), "--write"], cwd=ROOT, check=True
)
subprocess.run([sys.executable, str(ROOT / "comparison.py")], cwd=ROOT, check=True)
