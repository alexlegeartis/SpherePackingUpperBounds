"""Freeze successful radial exports, without registering unverified bounds."""

import argparse, json, hashlib, subprocess
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = "6e7990802b8e045a78e3f54c99669d5fe9b5c437"
baseline = json.loads(
    subprocess.check_output(
        ["git", "show", BASE + ":RESULTS.json"], cwd=ROOT, text=True
    )
)
WRAPPER = '''"""Replay every gate of this frozen certificate; assertions are required."""
import hashlib,json,subprocess,sys
from fractions import Fraction
from pathlib import Path
if not __debug__:raise RuntimeError('Assertions must remain enabled')
D=Path(__file__).resolve().parent
cfg=json.loads((D/'config.json').read_text())
for name,want in json.loads((D/'source_hashes.json').read_text()).items():
    assert hashlib.sha256((D/name).read_bytes()).hexdigest()==want,name
row=next(r for r in json.loads((D.parents[2]/'RESULTS.json').read_text()) if r['dimension']==cfg['dimension'])
if row['package']==cfg['package']:
    assert Fraction(row['strict_density_upper_bound'])==Fraction(cfg['strict_density_upper_bound'])
else:
    assert cfg['package'] in row.get('earlier_repository_packages',[])
    assert Fraction(row['strict_density_upper_bound'])<=Fraction(cfg['strict_density_upper_bound'])
executed=[]
for script in ['primary_angular.py','primary_radial.py','independent_angular.py','independent_radial.py']:
    result=subprocess.run([sys.executable,str(D.parents[1]/'profile_enrichment_20261114'/script),str(D)],capture_output=True,text=True)
    (D/(script.removesuffix('.py')+'.log')).write_text(result.stdout+result.stderr)
    print(result.stdout,flush=True)
    assert result.returncode==0,(script,result.stderr)
    executed.append({'script':script,'returncode':result.returncode})
out={'pass':True,'dimension':cfg['dimension'],'strict_density_upper_bound':cfg['strict_density_upper_bound'],
     'complete_primary_and_independent_replay':True,'executed':executed,
     'source_manifest_sha256':hashlib.sha256((D/'source_hashes.json').read_bytes()).hexdigest(),
     'execution':'Registered wrapper, including registry/hash gates and all four mathematical proof programs',
     'external_neighbor_counts_required':False,'arithmetic_backend':'FLINT/Arb shared; independent proof implementations'}
(D/'full_replay_result.json').write_text(json.dumps(out,indent=2)+'\\n')
print('FULL CERTIFICATE PASS',cfg['dimension'],cfg['strict_density_upper_bound'],flush=True)
'''
ap = argparse.ArgumentParser()
ap.add_argument("candidates", nargs="+")
args = ap.parse_args()
for candidate in args.candidates:
    source = Path(candidate).resolve()
    if source.is_dir():
        source = source / "radial_certificate.json"
    elif source.name == "candidate.json":
        source = source.with_name("radial_certificate.json")
    r = json.loads(source.read_text())
    assert r["radial_checks_pass"]
    n = r["dimension"]
    prior = next(z for z in baseline if z["dimension"] == n)
    assert Fraction(r["strict_density_upper_bound"]) < Fraction(
        prior["strict_density_upper_bound"]
    )
    D = ROOT / f"verifications/dim{n}/refinement_larger20261115"
    D.mkdir(parents=True, exist_ok=True)
    angular = json.loads((ROOT / r["angular_file"]).read_text())
    for key in ["gegenbauer_coefficients", "polynomial_coefficients"]:
        while Fraction(angular[key][-1]) == 0:
            angular[key].pop()
    angular["degree"] = len(angular["gegenbauer_coefficients"]) - 1
    angular["profile_a_exact"] = r["a_exact"]
    angular["profile_b_exact"] = r["b_exact"]
    angular.setdefault("negative_cap_delta", "1/100000")
    angular.setdefault("cutoff", "5001/10000")
    angular.setdefault("tail_u", 8)
    assert Fraction(angular["B_exact"]) == Fraction(r["B_exact"])
    for key in [
        "certified",
        "radial_checks_pass",
        "improves_registered",
        "scope",
        "pair_intervals",
        "pair_depth",
        "pair_seconds",
    ]:
        r.pop(key, None)
    r["data_status"] = (
        "Frozen rational witness; run verify.py for the density conclusion"
    )
    r["research_source"] = str(source.relative_to(ROOT))
    cfg = {
        "dimension": n,
        "package": str(D.relative_to(ROOT)),
        "strict_density_upper_bound": r["strict_density_upper_bound"],
        "previous_repository_strict_upper": prior["strict_density_upper_bound"],
        "previous_repository_package": prior["package"],
        "base_commit": BASE,
        "minimum_squared_distance": r["s_exact"],
        "profile_a_exact": r["a_exact"],
        "profile_b_exact": r["b_exact"],
        "radial_correction_degree": r["correction_degree"],
        "external_neighbor_count_inputs": False,
    }
    for name, data in [
        ("certificate.json", r),
        ("angular.json", angular),
        ("config.json", cfg),
    ]:
        (D / name).write_text(json.dumps(data, indent=2) + "\n")
    (D / "verify.py").write_text(WRAPPER)
    (D / "README.md").write_text(
        f"""# {n}D larger-profile refinement, 15 November 2026

Strict density bound: **{r['strict_density_upper_bound']}**.
Preceding main: {prior['strict_density_upper_bound']} at `{BASE}`.

The radial profile is `h(v)=v^(-a) exp(-b(v-1))`, with
`a={r['a_exact']}`, `b={r['b_exact']}`, and minimum squared distance
`s={r['s_exact']}`. The physical and Fourier polynomials have degree
{r['degree']}; the discovery correction space has degree {r['correction_degree']}.
The package uses no external neighbor counts.

See [the shared proof](../../profile_enrichment_20261114/PROOF.md).
Run `python verify.py` with python-flint installed and assertions enabled.
The wrapper checks the registry and source hashes, then executes all four
primary and independent angular/radial checks. Both paths share FLINT/Arb.
Discovery solver success alone does not certify a bound.
"""
    )
    names = ["certificate.json", "angular.json", "config.json", "verify.py"]
    names += [
        "../../profile_enrichment_20261114/" + p.name
        for p in sorted(
            (ROOT / "verifications/profile_enrichment_20261114").glob("*.py")
        )
    ]
    hashes = {
        name: hashlib.sha256((D / name).read_bytes()).hexdigest() for name in names
    }
    (D / "source_hashes.json").write_text(json.dumps(hashes, indent=2) + "\n")
    print(n, r["strict_density_upper_bound"], str(D.relative_to(ROOT)))
